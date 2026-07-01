import Base: eltype

# PharmacophoreGMM

struct PharmacophoreGMM{N,T<:Real,K,G<:SimpleMolGraph} <: AbstractIsotropicMultiGMM{N,T,K}
    gmms::Dict{K, IsotropicGMM{N,T}}
    graph::G
    σfun::Function
    ϕfun::Function
    feature_maps::Dict{K, Vector{Vector{Int}}}
end

eltype(::Type{PharmacophoreGMM{N,T,K,G}}) where {N,T,K,G} = Pair{K, IsotropicGMM{N,T}}

"""
    PharmacophoreGMM(mol; σfun=vdw_volume_sigma, ϕfun=a->one(typeof(vdw_radius(a))), feature_maps=…)

Build a `PharmacophoreGMM` from a molecule or subgraph `mol`, with one `IsotropicGMM`
per molecular-feature family named in `feature_maps`.

`feature_maps` maps each feature family (a `Symbol`) to a list of node-index sets; each
set is collapsed into a single `IsotropicGaussian` feature by `atoms_to_feature`. By
default it places one single-atom `:Volume` feature on every heavy atom of `mol`.
`feature_maps` derives pharmacophore families (donors, acceptors, rings, …) from feature
definitions.

`ϕfun(atom)` gives an atom's amplitude and `σfun(atom, ϕ)` its width; see
`atoms_to_feature` for how they combine over multi-atom feature sets.
"""
function PharmacophoreGMM(mol::SDFMolGraph;
                          σfun = vdw_volume_sigma,
                          ϕfun = a -> one(typeof(vdw_radius(a))),
                          feature_maps::Dict{K,Vector{Vector{Int}}} = Dict{Symbol,Vector{Vector{Int}}}(:Volume => [[i] for i in heavy_atom_idxs(mol)])) where K
    N = length(props(mol,1).coords)
    T = eltype(props(mol,1).coords)
    # add a GMM for each type of feature
    gmms = Dict{Symbol,IsotropicGMM{N,T}}()
    for (feature, nodesets) in feature_maps
        feats = IsotropicGaussian{N,T}[]
        for set in nodesets
            push!(feats, atoms_to_feature(mol, set; σfun=σfun, ϕfun=ϕfun))
        end
        push!(gmms, Pair(feature, IsotropicGMM(feats)))
    end
    return PharmacophoreGMM(gmms, mol, σfun, ϕfun, feature_maps)
end

function  Base.:*(R::AbstractMatrix{W}, x::PharmacophoreGMM{N,V,K,G}) where {N,V,K,G,W}
    numtype = promote_type(V, W)
    gmmdict = Dict{K, IsotropicGMM{N,numtype}}()
    for (key, gmm) in x.gmms
        push!(gmmdict, key=>R*gmm)
    end
    return PharmacophoreGMM{N,numtype,K,G}(gmmdict, x.graph, x.σfun, x.ϕfun, x.feature_maps)
end

function  Base.:+(x::PharmacophoreGMM{N,V,K,G}, T::AbstractVector{W}) where {N,V,K,G,W}
    numtype = promote_type(V, W)
    gmmdict = Dict{K, IsotropicGMM{N,numtype}}()
    for  (key, gmm) in x.gmms
        push!(gmmdict, key=>gmm+T)
    end
    return PharmacophoreGMM{N,numtype,K,G}(gmmdict, x.graph, x.σfun, x.ϕfun, x.feature_maps)
end

Base.:-(x::PharmacophoreGMM, T::AbstractVector) = x + (-T)

function transform(pgmm::PharmacophoreGMM{N,T,K,G}, tform) where {N,T,K,G}
    newgmms = Dict{K, IsotropicGMM{N,T}}()
    for (key, gmm) in pgmm.gmms
        push!(newgmms, key => tform(gmm))
    end
    return PharmacophoreGMM{N,T,K,G}(newgmms, transformed(tform, pgmm.graph), pgmm.σfun, pgmm.ϕfun, pgmm.feature_maps)
end

# `public` is a soft keyword only on Julia 1.11+; guard so the module still
# loads on the 1.10 LTS, where the declaration is simply absent.
@static if VERSION >= v"1.11"
    eval(Meta.parse("public transform"))
end

# descriptive display

Base.show(io::IO, pgmm::PharmacophoreGMM) = print(io,
    summary(pgmm),
    " with $(sum(length(gmm.second) for gmm in pgmm)) Gaussians in $(length(pgmm)) GMMs"
)

Base.show(io::IO, ::MIME"text/plain", pgmm::PharmacophoreGMM) = print(io,
    summary(pgmm),
    " from molecule with formula $(molecular_formula(pgmm.graph))",
    " with $(sum(length(gmm.second) for gmm in pgmm)) Gaussians in $(length(pgmm)) GMMs with labels:\n",
    "$([label for (label, gmm) in pgmm.gmms])"
)
