import Base: eltype

# PharmacophoreGMM

"""
    PharmacophoreGMM

A labeled isotropic Gaussian mixture built from a molecule. Every Gaussian in
`gaussians` carries a parallel feature label in `labels`, so overlap between two
`PharmacophoreGMM`s is restricted to Gaussian pairs whose labels interact (by
default, only equal labels; see [`overlap`](@ref) for the `interactions`
keyword). The two vectors have equal length.

The molecular `graph` is kept alongside the labeled vectors so the molecule can
still be drawn and its formula reported. `σfun`, `ϕfun`, and `feature_maps`
record how the Gaussians were built.
"""
struct PharmacophoreGMM{N,T<:Real,K,G<:SimpleMolGraph} <: AbstractLabeledIsotropicGMM{N,T,K}
    gaussians::Vector{IsotropicGaussian{N,T}}
    labels::Vector{K}
    graph::G
    σfun::Function
    ϕfun::Function
    feature_maps::Dict{K, Vector{Vector{Int}}}
    function PharmacophoreGMM{N,T,K,G}(gaussians, labels, graph, σfun, ϕfun, feature_maps) where {N,T,K,G}
        length(gaussians) == length(labels) ||
            throw(DimensionMismatch("number of Gaussians ($(length(gaussians))) must match number of labels ($(length(labels)))"))
        return new{N,T,K,G}(gaussians, labels, graph, σfun, ϕfun, feature_maps)
    end
end

function PharmacophoreGMM(gaussians::AbstractVector{IsotropicGaussian{N,T}}, labels::AbstractVector{K},
                          graph::G, σfun, ϕfun, feature_maps) where {N,T,K,G<:SimpleMolGraph}
    return PharmacophoreGMM{N,T,K,G}(gaussians, labels, graph, σfun, ϕfun, feature_maps)
end

eltype(::Type{PharmacophoreGMM{N,T,K,G}}) where {N,T,K,G} = IsotropicGaussian{N,T}

"""
    feature_labels(pgmm::PharmacophoreGMM)

Return the distinct feature labels present in `pgmm`, in first-appearance order.
"""
feature_labels(pgmm::PharmacophoreGMM) = unique(pgmm.labels)

"""
    PharmacophoreGMM(mol; σfun=vdw_volume_sigma, ϕfun=a->one(typeof(vdw_radius(a))), feature_maps=…)

Build a `PharmacophoreGMM` from a molecule or subgraph `mol`: one labeled
`IsotropicGaussian` per molecular-feature set named in `feature_maps`.

`feature_maps` maps each feature family (a `Symbol`) to a list of node-index sets; each
set is collapsed into a single `IsotropicGaussian` feature by `atoms_to_feature` and
labeled with its family. By default it places one single-atom `:Volume` feature on every
heavy atom of `mol`. `feature_maps` derives pharmacophore families (donors, acceptors,
rings, …) from feature definitions.

`ϕfun(atom)` gives an atom's amplitude and `σfun(atom, ϕ)` its width; see
`atoms_to_feature` for how they combine over multi-atom feature sets.

# Example

```jldoctest
julia> mol = sdftomol(joinpath(pkgdir(MolecularGaussians), "assets", "data", "E1050_3d.sdf"));

julia> PharmacophoreGMM(mol)
PharmacophoreGMM{3, Float64, Symbol, SDFMolGraph} from molecule with formula C18H24O8S2 with 28 Gaussians labeled:
[:Volume]
```
"""
function PharmacophoreGMM(mol::SDFMolGraph;
                          σfun = vdw_volume_sigma,
                          ϕfun = a -> one(typeof(vdw_radius(a))),
                          feature_maps::Dict{K,Vector{Vector{Int}}} = Dict{Symbol,Vector{Vector{Int}}}(:Volume => [[i] for i in heavy_atom_idxs(mol)])) where K
    N = length(props(mol,1).coords)
    T = eltype(props(mol,1).coords)
    gaussians = IsotropicGaussian{N,T}[]
    labels = K[]
    for (feature, nodesets) in feature_maps
        for set in nodesets
            push!(gaussians, atoms_to_feature(mol, set; σfun=σfun, ϕfun=ϕfun))
            push!(labels, feature)
        end
    end
    return PharmacophoreGMM(gaussians, labels, mol, σfun, ϕfun, feature_maps)
end

# GaussianMixtureAlignment's generic `*`/`+` on AbstractLabeledIsotropicGMM
# return a bare LabeledIsotropicGMM, dropping the graph and feature metadata.
# These methods keep the PharmacophoreGMM type and carry those fields through.
function  Base.:*(R::AbstractMatrix{W}, x::PharmacophoreGMM{N,V,K,G}) where {N,V,K,G,W}
    numtype = promote_type(V, W)
    gaussians = IsotropicGaussian{N,numtype}[R*g for g in x.gaussians]
    return PharmacophoreGMM{N,numtype,K,G}(gaussians, x.labels, x.graph, x.σfun, x.ϕfun, x.feature_maps)
end

function  Base.:+(x::PharmacophoreGMM{N,V,K,G}, T::AbstractVector{W}) where {N,V,K,G,W}
    numtype = promote_type(V, W)
    gaussians = IsotropicGaussian{N,numtype}[g+T for g in x.gaussians]
    return PharmacophoreGMM{N,numtype,K,G}(gaussians, x.labels, x.graph, x.σfun, x.ϕfun, x.feature_maps)
end

Base.:-(x::PharmacophoreGMM, T::AbstractVector) = x + (-T)

"""
    transform(pgmm::PharmacophoreGMM, tform) -> PharmacophoreGMM

Apply the transformation `tform` to every Gaussian of `pgmm` and to its
underlying molecular graph, returning a new `PharmacophoreGMM`. `tform` is
called as `tform(g)` on each `IsotropicGaussian` and must map a set of 3-D
points rigidly (e.g. an `AffineMap` from CoordinateTransformations).
"""
function transform(pgmm::PharmacophoreGMM, tform)
    gaussians = [tform(g) for g in pgmm.gaussians]
    return PharmacophoreGMM(gaussians, pgmm.labels, transformed(tform, pgmm.graph), pgmm.σfun, pgmm.ϕfun, pgmm.feature_maps)
end

# `public` is a soft keyword only on Julia 1.11+; guard so the module still
# loads on the 1.10 LTS, where the declaration is simply absent.
@static if VERSION >= v"1.11"
    eval(Meta.parse("public transform"))
end

# descriptive display

Base.show(io::IO, pgmm::PharmacophoreGMM) = print(io,
    summary(pgmm),
    " with $(length(pgmm.gaussians)) Gaussians labeled $(feature_labels(pgmm))"
)

Base.show(io::IO, ::MIME"text/plain", pgmm::PharmacophoreGMM) = print(io,
    summary(pgmm),
    " from molecule with formula $(molecular_formula(pgmm.graph))",
    " with $(length(pgmm.gaussians)) Gaussians labeled:\n",
    "$(feature_labels(pgmm))"
)
