import Base: eltype

# AtomGaussian

struct AtomGaussian{N,T<:Real} <: AbstractIsotropicGaussian{N,T}
    μ::SVector{N,T}
    σ::T
    ϕ::T
    node::SDFileAtom
    nodeidx::Int
end

"""
"""
function AtomGaussian(atom::SDFileAtom, nodeidx, σ, ϕ)
    N = length(atom.coords)
    T = eltype(atom.coords)
    return AtomGaussian{N,T}(SVector{N,T}(atom.coords), T(σ), T(ϕ), atom, nodeidx)
end

AtomGaussian(mol::Union{UndirectedGraph, SubgraphView}, nodeidx, σdict, ϕdict
    ) = AtomGaussian(nodeattr(mol,nodeidx), nodeidx, σdict[nodeidx], ϕdict[nodeidx])

function AtomGaussian(μ::AbstractArray, σ::Real, ϕ::Real, node=SDFileAtom(), nodeidx=0)
    t = promote_type(eltype(μ), typeof(σ), typeof(ϕ))
    return AtomGaussian{length(μ),t}(SVector{length(μ),t}(μ), t(σ), t(ϕ), node, nodeidx)
end

convert(::Type{AtomGaussian{N,T}}, g::AtomGaussian) where {N,T} = AtomGaussian(SVector{N,T}(g.μ), T(g.σ), T(g.ϕ), g.node, g.nodeidx)

function Base.:*(R::AbstractMatrix{W}, x::AtomGaussian{N,V}) where {N,V,W}
    numtype = promote_type(V, W)
    return AtomGaussian{N,numtype}(R*x.μ, x.σ, x.ϕ, x.node, x.nodeidx)
end

function Base.:+(x::AtomGaussian{N,V}, T::AbstractVector{W}) where {N,V,W}
    numtype = promote_type(V, W)
    return AtomGaussian{N,numtype}(x.μ.+T, x.σ, x.ϕ, x.node, x.nodeidx)
end

Base.:-(x::AtomGaussian, T::AbstractVector,) = x + (-T)


# MolGMM

struct MolGMM{N,T<:Real,G<:Union{UndirectedGraph,SubgraphView}} <: AbstractIsotropicGMM{N,T}
    gaussians::Vector{AtomGaussian{N,T}}
    graph::G
    nodes::Set{Int}
    σfun::Function
    ϕfun::Function
end

eltype(::Type{MolGMM{N,T}}) where {N,T} = AtomGaussian{N,T}

function gauss_from_atom(mol::UndirectedGraph, index::Int, σfun, ϕfun)
    atom = nodeattr(mol, index)
    return AtomGaussian(atom, index, σfun(atom), ϕfun(atom))
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
                nodes=nodeset(mol);
                σfun = vdw_volume_sigma,
                ϕfun = rocs_volume_amplitude)
    N = length(nodeattrs(mol)[1].coords)
    T = eltype(nodeattrs(mol)[1].coords)
    atoms = [gauss_from_atom(mol, i, σfun, ϕfun) for i in nodes]
    # remove atoms with ϕ = 0
    filter!(atom->atom.ϕ≠0,atoms)
    return MolGMM{N,T,typeof(mol)}(atoms, mol, nodes, σfun, ϕfun)
end

MolGMM(submol::SubgraphView; kwargs...) = MolGMM(submol.graph, nodeset(submol); kwargs...)

function Base.:*(R::AbstractMatrix{W}, x::MolGMM{N,V,G}) where {N,V,G,W}
    numtype = promote_type(V, W)
    return MolGMM{N,numtype,G}([R*g for g in x.gaussians], x.graph, x.nodes, x.σfun, x.ϕfun)
end

function Base.:+(x::MolGMM{N,V,G}, T::AbstractVector{W}) where {N,V,G,W}
    numtype = promote_type(V, W)
    return MolGMM{N,numtype,G}([g+T for g in x.gaussians], x.graph, x.nodes, x.σfun, x.ϕfun)
end

Base.:-(x::MolGMM, T::AbstractVector,) = x + (-T)


# PharmacophoreGMM

struct PharmacophoreGMM{N,T<:Real,K,G<:Union{UndirectedGraph,SubgraphView}} <: AbstractIsotropicMultiGMM{N,T,K}
    gmms::Dict{K, IsotropicGMM{N,T}}
    graph::G
    nodes::Set{Int}
    σfun::Function
    ϕfun::Function
    features::Vector{K}
end

eltype(::Type{PharmacophoreGMM{N,T,K,G}}) where {N,T,K,G} = Pair{K, IsotropicGMM{N,T,G}}

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
                          features = pubchem_features)
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
                    push!(feats, atoms_to_feature(mol, set; σfun=σfun, ϕfun=ϕfun))
                end
            end
            push!(gmms, Pair(key, IsotropicGMM(feats)))
        end
    end
    gmms = gmms
    return PharmacophoreGMM(gmms, mol, nodes, σfun, ϕfun, features)
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
