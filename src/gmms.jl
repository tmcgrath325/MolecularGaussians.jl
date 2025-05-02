import Base: eltype

# PharmacophoreGMM

struct PharmacophoreGMM{N,T<:Real,K,G<:SimpleMolGraph} <: GaussianMixtureAlignment.AbstractLabeledIsotropicGMM{N,T,K}
    gaussians::Vector{LabeledIsotropicGaussian{N,T,K}}
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
                          ϕfun = a -> one(typeof(MolecularGraph.atom_radius(a)));
                          featuremaps::Dict{K,Vector{Vector{Int}}} = Dict{Symbol,Vector{Vector{Int}}}(:Volume => [[i] for i in heavy_atom_idxs(mol)])) where K
    N = length(props(mol,1).coords)
    T = eltype(props(mol,1).coords)
    # add a GMM for each type of feature
    gaussians = Vector{LabeledIsotropicGaussian{N,T,Symbol}}()
    for (feature, nodesets) in featuremaps
        for set in nodesets
            push!(gaussians, atoms_to_feature(mol, set, feature; σfun=σfun, ϕfun=ϕfun))
        end
    end
    return PharmacophoreGMM(gaussians, mol, σfun, ϕfun, featuremaps)
end

function  Base.:*(R::AbstractMatrix{W}, x::PharmacophoreGMM{N,V,K,G}) where {N,V,K,G,W}
    numtype = promote_type(V, W)
    return PharmacophoreGMM{N,numtype,K,G}([R*g for g in x.gaussians], x.graph, x.σfun, x.ϕfun, x.featuremaps)
end

function  Base.:+(x::PharmacophoreGMM{N,V,K,G}, T::AbstractVector{W}) where {N,V,K,G,W}
    numtype = promote_type(V, W)
    return PharmacophoreGMM{N,numtype,K,G}([g+T for g in x.gaussians], x.graph, x.σfun, x.ϕfun, x.featuremaps)
end

Base.:-(x::PharmacophoreGMM, T::AbstractVector) = x + (-T)

# function transform(pgmm::PharmacophoreGMM{N,T,K,G}, tform::AffineMap) where {N,T,K,G}
#     return PharmacophoreGMM{N,T,K,G}([tform(g) for g in pgmm.gaussians], tform(pgmm.graph), pgmm.σfun, pgmm.ϕfun, pgmm.featuremaps)
# end

feature_labels(x::PharmacophoreGMM) = unique([g.label for g in x.gaussians])

# descriptive display

Base.show(io::IO, pgmm::PharmacophoreGMM) = println(io,
    summary(pgmm),
    " from molecule with formula $(molecular_formula(pgmm.graph))",
    " with $(length(pgmm)) Gaussians with labels:\n",
    "$([label for label in feature_labels(pgmm)])"
)
