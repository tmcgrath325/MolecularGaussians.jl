import Base: eltype

# AtomGaussian

struct AtomGaussian{N,T<:Real,A} <: AbstractIsotropicGaussian{N,T}
    μ::SVector{N,T}
    σ::T
    ϕ::T
    node::A
    nodeidx::Int
end

"""
"""
function AtomGaussian(atom::SDFAtom, nodeidx, σ, ϕ)
    N = length(atom.coords)
    T = eltype(atom.coords)
    return AtomGaussian{N,T,SDFAtom}(SVector{N,T}(atom.coords), T(σ), T(ϕ), atom, nodeidx)
end

AtomGaussian(mol::SimpleMolGraph, nodeidx, σdict, ϕdict
    ) = AtomGaussian(props(mol,nodeidx), nodeidx, σdict[nodeidx], ϕdict[nodeidx])

function AtomGaussian(μ::AbstractArray, σ::Real, ϕ::Real, node=SDFAtom(), nodeidx=0)
    t = promote_type(eltype(μ), typeof(σ), typeof(ϕ))
    return AtomGaussian{length(μ),t}(SVector{length(μ),t}(μ), t(σ), t(ϕ), node, nodeidx)
end

convert(::Type{AtomGaussian{N,T}}, g::AtomGaussian) where {N,T} = AtomGaussian(SVector{N,T}(g.μ), T(g.σ), T(g.ϕ), g.node, g.nodeidx)

function Base.:*(R::AbstractMatrix{W}, x::AtomGaussian{N,V,A}) where {N,V,A,W}
    numtype = promote_type(V, W)
    return AtomGaussian{N,numtype,A}(R*x.μ, x.σ, x.ϕ, x.node, x.nodeidx)
end

function Base.:+(x::AtomGaussian{N,V,A}, T::AbstractVector{W}) where {N,V,A,W}
    numtype = promote_type(V, W)
    return AtomGaussian{N,numtype,A}(x.μ.+T, x.σ, x.ϕ, x.node, x.nodeidx)
end

Base.:-(x::AtomGaussian, T::AbstractVector,) = x + (-T)


# MolGMM

struct MolGMM{N,T<:Real,A,G<:SimpleMolGraph} <: AbstractIsotropicGMM{N,T}
    gaussians::Vector{AtomGaussian{N,T,A}}
    graph::G
    σfun::Function
    ϕfun::Function
end

eltype(::Type{MolGMM{N,T,A}}) where {N,T,A} = AtomGaussian{N,T,A}

function gauss_from_atom(mol::SimpleMolGraph, index::Int, σfun, ϕfun)
    atom = props(mol, index)
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
function MolGMM(mol::SimpleMolGraph,
                σfun = vdw_volume_sigma,
                ϕfun = rocs_volume_amplitude,
                nodes = vertices(mol))
    N = length(props(mol,1).coords)
    T = eltype(props(mol,1).coords)
    A = eltype(values(vprops(mol)))
    atoms = [gauss_from_atom(mol, i, σfun, ϕfun) for i in nodes]
    # remove atoms with ϕ = 0
    filter!(atom->atom.ϕ≠0,atoms)
    return MolGMM{N,T,A,typeof(mol)}(atoms, mol, σfun, ϕfun)
end

function Base.:*(R::AbstractMatrix{W}, x::MolGMM{N,V,A,G}) where {N,V,A,G,W}
    numtype = promote_type(V, W)
    return MolGMM{N,numtype,A,G}([R*g for g in x.gaussians], x.graph, x.σfun, x.ϕfun)
end

function Base.:+(x::MolGMM{N,V,A,G}, T::AbstractVector{W}) where {N,V,A,G,W}
    numtype = promote_type(V, W)
    return MolGMM{N,numtype,A,G}([g+T for g in x.gaussians], x.graph, x.σfun, x.ϕfun)
end

Base.:-(x::MolGMM, T::AbstractVector,) = x + (-T)


# PharmacophoreGMM

struct PharmacophoreGMM{N,T<:Real,K,G<:SimpleMolGraph} <: AbstractIsotropicMultiGMM{N,T,K}
    gmms::Dict{K, IsotropicGMM{N,T}}
    graph::G
    σfun::Function
    ϕfun::Function
    featuremaps::Dict{K, Vector{Vector{Int}}}
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
function PharmacophoreGMM(mol::SDFMolGraph,
                          σfun = vdw_volume_sigma,
                          ϕfun = rocs_volume_amplitude;
                          featuremaps::Dict{K,Vector{Vector{Int}}} = Dict{Symbol,Vector{Vector{Int}}}(:Volume => [[i] for i in heavy_atom_idxs(mol)])) where K
    N = length(props(mol,1).coords)
    T = eltype(props(mol,1).coords)
    # add a GMM for each type of feature
    gmms = Dict{Symbol,IsotropicGMM{N,T}}()
    for (feature, nodesets) in featuremaps
        feats = IsotropicGaussian{N,T}[]
        for set in nodesets
            push!(feats, atoms_to_feature(mol, set; σfun=σfun, ϕfun=ϕfun))
        end
        push!(gmms, Pair(feature, IsotropicGMM(feats)))
    end
    return PharmacophoreGMM(gmms, mol, σfun, ϕfun, featuremaps)
end

function  Base.:*(R::AbstractMatrix{W}, x::PharmacophoreGMM{N,V,K,G}) where {N,V,K,G,W}
    numtype = promote_type(V, W)
    gmmdict = Dict{K, IsotropicGMM{N,numtype}}()
    for (key, gmm) in x.gmms
        push!(gmmdict, key=>R*gmm)
    end
    return PharmacophoreGMM{N,numtype,K,G}(gmmdict, R*x.graph, x.nodes, x.σfun, x.ϕfun, x.features)
end

function  Base.:+(x::PharmacophoreGMM{N,V,K,G}, T::AbstractVector{W}) where {N,V,K,G,W}
    numtype = promote_type(V, W)
    gmmdict = Dict{K, IsotropicGMM{N,numtype}}()
    for  (key, gmm) in x.gmms
        push!(gmmdict, key=>gmm+T)
    end
    return PharmacophoreGMM{N,numtype,K,G}(gmmdict, x.graph+T, x.nodes, x.σfun, x.ϕfun, x.features)
end

Base.:-(x::PharmacophoreGMM, T::AbstractVector) = x + (-T)

# descriptive display

Base.show(io::IO, molgmm::MolGMM) = println(io,
    "MolGMM from molecule with formula $(molecular_formula(molgmm.graph))",
    " with $(length(molgmm)) Gaussians."
)

Base.show(io::IO, pgmm::PharmacophoreGMM) = println(io,
    summary(pgmm),
    " from molecule with formula $(molecular_formula(pgmm.graph))",
    " with $(sum(length(gmm.second) for gmm in pgmm)) Gaussians in $(length(pgmm)) GMMs with labels:\n",
    "$([label for (label, gmm) in pgmm.gmms])"
)
