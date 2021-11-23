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

eltype(::Type{MolGMM{N,T}}) where {N,T} = AtomGaussian{N,T}

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
                nodes=nodeset(mol);
                σfun = vdwvolume_sigma,
                ϕfun = ones)
    N = length(nodeattrs(mol)[1].coords)
    T = eltype(nodeattrs(mol)[1].coords)
    σdict, ϕdict = Dict{Int, T}(σfun(mol)), Dict{Int, T}(ϕfun(mol))
    atoms = [AtomGaussian(nodeattr(mol, i), i, σdict[i], ϕdict[i]) for i in nodes]
    # remove atoms with ϕ = 0
    filter!(atom->atom.ϕ≠0,atoms)
    return MolGMM{N,T}(atoms, mol, nodes, σfun, ϕfun)
end

MolGMM(submol::SubgraphView; kwargs...) = MolGMM(submol.graph, nodeset(submol); kwargs...)

struct PharmacophoreGMM{N,T<:Real,K} <: AbstractIsotropicMultiGMM{N,T,K}
    gmms::Dict{K, IsotropicGMM{N,T}}
    graph::Union{UndirectedGraph,SubgraphView}
    nodes::Set{Int}
    σfun::Function
    ϕfun::Function
    features::Vector{K}
    directional::Bool
end

eltype(::Type{PharmacophoreGMM{N,T,K}}) where {N,T,K} = Pair{K, IsotropicGMM{N,T}}

"""
    model = PharmacophoreGMM(mol, σfun=vdwvolume_sigma, ϕfun=ones, nodes=nodeset(mol), features=pubchem_features)

Creates a set Gaussian mixture models from a molecule or subgraph `mol`, with each model corresponding to a
particular type of molecular feature (e.g. ring structures)

Optionally, functions `σfun` and `ϕfun` can be provided, which take `mol` as input and return dictionaries
mapping node indicies to variances `σ` and scaling coefficients `ϕ`, respectively.

If `nodes` is provided, the Gaussian mixture models will be constructed only from atoms corresponding to the node
indexes of the molecule's graph.

`features` specifies the keys of allowed molecular features. Features with other keys will be ignored.

`directional` specifies whether or not geometric constraints for ring structures and hydrogen bond donors will be included.
"""
function PharmacophoreGMM(mol::UndirectedGraph,
                          nodes = nodeset(mol);
                          σfun = vdwvolume_sigma,
                          ϕfun = ones,
                          features = pubchem_features,
                          directional = true)
    dim = length(nodeattrs(mol)[1].coords)
    valtype = eltype(nodeattrs(mol)[1].coords)
    # add a GMM for each type of feature
    gmms = Dict{Symbol,IsotropicGMM{dim,valtype}}()
    for (key, nodesets) in pharmfeatures!(mol)
        if key in features
            feats = IsotropicGaussian{dim,valtype}[]
            for set in nodesets
                # check if all atoms in the feature are represented by allowed nodes
                if set ∩ nodes == set
                    # if it is a ring feature, add normal vectors
                    if (key == :rings || key == :aromaticrings) && directional
                        center = sum([nodeattr(mol,idx).coords for idx in set])/length(set)
                        coordmat = fill(zero(valtype), 3, length(set))
                        for (i,idx) in enumerate(set)
                            coordmat[:,i] = nodeattr(mol,idx).coords .- center
                        end
                        n = GaussianMixtureAlignment.planefit(coordmat)[1]
                        dirs = [n, -n]
                    # if it is a hydrogen bond donor, point outward from hydrogens
                    elseif key == :donor && directional
                        if length(set)==1
                            dirs = Vector{valtype}[]
                            idx = collect(set)[1]
                            for (edge,neigh) in neighbors(mol,idx)
                                if nodeattr(mol,neigh).symbol == :H
                                    push!(dirs, nodeattr(mol,neigh).coords - nodeattr(mol,idx).coords)
                                end
                            end
                        else
                            dirs = nothing
                        end
                    # if it is a hydrogen bond acceptor, point along the lone pair orbitals
                    elseif key == :acceptor && directional
                        # the acceptor should be a single atom
                        if length(set) != 1
                            dirs = nothing
                        else
                            idx = collect(set)[1]
                            neighs = collect(neighbors(mol,idx))
                            if hybridization(mol)[idx] == :sp3      # tetrahedral
                                if length(neighs) == 3
                                    # if there are 3 bonds
                                    bonddir1 = nodeattr(mol,neighs[1].second).coords .- nodeattr(mol,idx).coords
                                    bonddir2 = nodeattr(mol,neighs[2].second).coords .- nodeattr(mol,idx).coords
                                    bonddir3 = nodeattr(mol,neighs[3].second).coords .- nodeattr(mol,idx).coords
                                    dirs = Vector{valtype}[-normalize(bonddir1.+bonddir2.+bonddir3)]
                                elseif length(neighs) == 2
                                    # if there are 2 bonds
                                    bonddir1 = nodeattr(mol,neighs[1].second).coords .- nodeattr(mol,idx).coords
                                    bonddir2 = nodeattr(mol,neighs[2].second).coords .- nodeattr(mol,idx).coords
                                    axis1 = normalize(bonddir1.+bonddir2)
                                    axis2 = normalize(cross(bonddir1,bonddir2))
                                    R = AngleAxis(π,axis2...) * AngleAxis(π/2,axis1...)
                                    dirs = Vector{valtype}[R*bonddir1, R*bonddir2]
                                else
                                    dirs = nothing
                                end
                            elseif hybridization(mol)[idx] == :sp2  # trigonal planar
                                # if there are two bonds, use them to find the lone pair's direction
                                if length(neighs) != 2
                                    dirs = nothing
                                else
                                    dirs = Vector{valtype}[]
                                    bonddir1 = nodeattr(mol,neighs[1].second).coords .- nodeattr(mol,idx).coords
                                    bonddir2 = nodeattr(mol,neighs[2].second).coords .- nodeattr(mol,idx).coords
                                    axis = normalize(cross(bonddir1, bonddir2))
                                    push!(dirs, AngleAxis(4π/3, axis...)*bonddir1)
                                end
                            elseif hybridization(mol)[idx] == :sp   # linear
                                # direction points along the axis of the atom's single bond
                                if length(neighs) != 1
                                    dirs = nothing
                                else
                                    dirs = Vector{valtype}[nodeattr(nodeattr(mol,idx) .- mol,neighs[1].second.coords)]
                                end
                            end
                        end
                    # otherwise there is no geometric constraint
                    else
                        dirs = nothing
                    end
                    push!(feats, atoms_to_feature(mol, set, dirs; σfun=σfun, ϕfun=ϕfun))
                end
            end
            push!(gmms, Pair(key, IsotropicGMM(feats)))
        end
    end
    gmms = gmms
    return PharmacophoreGMM(gmms, mol, nodes, σfun, ϕfun, features, directional)
end

PharmacophoreGMM(submol::SubgraphView; kwargs...) = PharmacophoreGMM(submol.graph, nodeset(submol); kwargs...)

# descriptive display

Base.show(io::IO, molgmm::MolGMM) = println(io,
    "MolGMM from molecule with formula $(typeof(molgmm.graph)<:GraphMol ? molecularformula(molgmm.graph) : molecularformula(molgmm.graph.graph))",
    " with $(length(molgmm)) Gaussians."
)

Base.show(io::IO, pgmm::PharmacophoreGMM) = println(io,
    summary(pgmm),
    " from molecule with formula $(typeof(pgmm.graph)<:GraphMol ? molecularformula(pgmm.graph) : molecularformula(pgmm.graph.graph))",
    " with $(sum(length(gmm.second) for gmm in pgmm)) Gaussians in $(length(pgmm)) GMMs with labels:\n",
    "$([label for (label, gmm) in pgmm.gmms])"
)
