import Base: eltype

# AtomGaussian

struct AtomGaussian{N,T<:Real} <: AbstractIsotropicGaussian{N,T}
    μ::SVector{N,T}
    σ::T
    ϕ::T
    dirs::Vector{SVector{N,T}}
    node::SDFileAtom
    nodeidx::Int
end

"""
"""
function AtomGaussian(atom::SDFileAtom, nodeidx, σ, ϕ)
    N = length(atom.coords)
    T = eltype(atom.coords)
    return AtomGaussian{N,T}(SVector{N,T}(atom.coords), T(σ), T(ϕ), SVector{N,T}[], atom, nodeidx)
end

AtomGaussian(mol::Union{UndirectedGraph, SubgraphView}, nodeidx, σdict, ϕdict
    ) = AtomGaussian(nodeattr(mol,nodeidx), nodeidx, σdict[nodeidx], ϕdict[nodeidx])

function AtomGaussian(μ::AbstractArray, σ::Real, ϕ::Real, dirs::AbstractArray=SVector{length(μ),eltype(μ)}[], node=SDFileAtom(), nodeidx=0)
    t = promote_type(eltype(μ), typeof(σ), typeof(ϕ), eltype(eltype(dirs)))
    return AtomGaussian{length(μ),t}(SVector{length(μ),t}(μ), t(σ), t(ϕ), SVector{length(μ),t}[SVector{length(μ),t}(dir/norm(dir)) for dir in dirs], node, nodeidx)
end

convert(::Type{AtomGaussian{N,T}}, g::AtomGaussian) where {N,T} = AtomGaussian(SVector{N,T}(g.μ), T(g.σ), T(g.ϕ), Vector{SVector{N,T}}(g.dirs), g.node, g.nodeidx)


# MolGMM

struct MolGMM{N,T<:Real} <: AbstractIsotropicGMM{N,T}
    gaussians::Vector{AtomGaussian{N,T}}
    graph::Union{UndirectedGraph,SubgraphView}
    nodes::Set{Int}
    σfun::Function
    ϕfun::Function
end

"""
    model = MolGMM(mol, σfun=ones, ϕfun=ones, nodes=nodeset(mol))

Creates a Gaussian mixture model from a molecule or subgraph `mol`, which consists of spherical Gaussian distributions 
with means `μ` equal to atom coordinates. 

Optionally, functions `σfun` and `ϕfun` can be provided, which take `mol` as input and return dictionaries mapping
mapping node indicies to variances `σ` and scaling coefficients `ϕ`, respectively.

If `nodes` is provided, the Gaussian mixture model will be constructed only from atoms corresponding to the node
indexes of the molecule's graph.
"""
function MolGMM(mol::UndirectedGraph,
                σfun = ones, 
                ϕfun = ones,
                nodes=nodeset(mol))
    N = length(nodeattrs(mol)[1].coords)
    T = eltype(nodeattrs(mol)[1].coords)
    σdict, ϕdict = Dict{Int, T}(σfun(mol)), Dict{Int, T}(ϕfun(mol))
    atoms = [AtomGaussian(nodeattr(mol, i), i, σdict[i], ϕdict[i]) for i in nodes]
    # remove atoms with ϕ = 0
    filter!(atom->atom.ϕ≠0,atoms)
    return MolGMM{N,T}(atoms, mol, nodes, σfun, ϕfun)
end

eltype(::MolGMM{N,T}) where {N,T} = AtomGaussian{N,T}

MolGMM(submol::SubgraphView, σfun=ones, ϕfun=ones) = MolGMM(submol.graph, σfun, ϕfun, nodeset(submol))

struct PharmacophoreGMM{N,T<:Real,K} <: AbstractIsotropicMultiGMM{N,T,K}
    gmms::Dict{K, IsotropicGMM{N,T}}
    directional::Bool
    graph::Union{UndirectedGraph,SubgraphView}
    nodes::Set{Int}
    σfun::Function
    ϕfun::Function
end


"""
    model = PharmacophoreGMM(mol, σfun=vdwvolume_sigma, ϕfun=ones, nodes=nodeset(mol), features=pubchem_features)

Creates a set Gaussian mixture models from a molecule or subgraph `mol`, with each model corresponding to a 
particular type of molecular feature (e.g. ring structures)

Optionally, functions `σfun` and `ϕfun` can be provided, which take `mol` as input and return dictionaries mapping
mapping node indicies to variances `σ` and scaling coefficients `ϕ`, respectively.

If `nodes` is provided, the Gaussian mixture models will be constructed only from atoms corresponding to the node
indexes of the molecule's graph.

`features` specifies the keys of allowed molecular features. Features with other keys will be ignored.

`directional` specifies whether or not geometric constraints for ring structures and hydrogen bond donors will be included.
"""
function PharmacophoreGMM(mol::UndirectedGraph,
                       σfun = vdwvolume_sigma, # vdwradii!,
                       ϕfun = ones,
                       nodes = nodeset(mol);
                       features = pubchem_features,
                       directional = true)
    dim = length(nodeattrs(mol)[1].coords)
    valtype = eltype(nodeattrs(mol)[1].coords)
    # add a GMM for each type of feature
    gmms = Dict{Symbol,IsotropicGMM{dim,valtype}}()
    for p in pharmfeatures!(mol)
        if p.first in features
            feats = IsotropicGaussian{dim,valtype}[]
            for set in p.second
                # check if all atoms in the feature are represented by allowed nodes
                if set ∩ nodes == set
                    # if it is a ring feature, add normal vectors
                    if (p.first == :rings || p.first == :aromaticrings) && directional
                        center = sum([nodeattr(mol,idx).coords for idx in set])/length(set)
                        coordmat = fill(zero(valtype), 3, length(set))
                        for (i,idx) in enumerate(set)
                            coordmat[:,i] = nodeattr(mol,idx).coords .- center
                        end
                        n = GaussianMixtureAlignment.planefit(coordmat)[1]
                        dirs = [n, -n]
                    # if it is a hydrogen bond donor, point outward from hydrogens
                    elseif p.first == :donor && directional
                        dirs = Vector{valtype}[]   
                        for idx in set
                            for neigh in neighbors(mol,idx)
                                if nodeattr(mol,neigh.second).symbol == :H
                                    push!(dirs, nodeattr(mol,neigh.second).coords - nodeattr(mol,idx).coords)
                                end
                            end
                        end
                    # if it is a hydrogen bond acceptor, point along the lone pair orbitals
                    elseif p.first == :acceptor && directional
                        # the acceptor should be a single atom
                        if length(set) != 1
                            dirs = nothing
                        else
                            idx = collect(set)[1]
                            if hybridization(mol)[idx] == :sp3      # tetrahedral
                                # TODO
                                # need at least 2 bonds, then find the other directions
                                dirs = nothing
                            elseif hybridization(mol)[idx] == :sp2  # planar
                                # TODO
                                # need to get the plane using hybridization of the atom's neighbors (?)
                                dirs = nothing
                            end
                        end
                    # otherwise there is no geometric constraint
                    else
                        dirs = nothing
                    end
                    push!(feats, atoms_to_feature(mol, set, σfun, ϕfun, dirs))
                end
            end
            push!(gmms, Pair(p.first, IsotropicGMM(feats)))
        end
    end
    gmms = gmms
    return PharmacophoreGMM(gmms, directional, mol, nodes, σfun, ϕfun)
end

PharmacophoreGMM(submol::SubgraphView, σfun=ones, ϕfun=ones, features=pubchem_features, directional=true) = 
    PharmacophoreGMM(submol.graph, σfun, ϕfun, nodeset(submol); features=features,directional=directional)

# descriptive display

Base.show(io::IO, molgmm::MolGMM) = println(io,
    "MolGMM from molecule with formula $(typeof(molgmm.graph)<:GraphMol ? molecularformula(molgmm.graph) : molecularformula(molgmm.graph.graph))",
    " with $(length(molgmm)) Gaussians."
)

Base.show(io::IO, pgmm::PharmacophoreGMM) = println(io,
    summary(pgmm),
    " from molecule with formula $(typeof(pgmm.graph)<:GraphMol ? molecularformula(pgmm.graph) : molecularformula(pgmm.graph.graph))",
    " with $(length(pgmm)) Gaussians in $(length(pgmm)) GMMs with labels:\n",
    "$([label for (label, gmm) in pgmm.gmms])"
)
