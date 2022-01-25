const pubchem_features = [:acceptor, :anion, :cation, :donor, :hydrophobe, :rings, :aromaticrings]

mutable struct FeatureGaussian{N,T<:Real,K} <: AbstractIsotropicGaussian{N,T}
    μ::SVector{N,T}
    σ::T
    ϕ::T
    dirs::Vector{SVector{N,T}}
    nodes::Vector{SDFileAtom}
    nodeidxs::Set{Int}
    type::K
end

function FeatureGaussian(μ::AbstractArray, σ::Real, ϕ::Real, dirs::AbstractArray=SVector{length(μ),eltype(μ)}[], nodes=SDFileAtom[], nodeidxs=[Set{Int}()], type=nothing)
    t = promote_type(eltype(μ), typeof(σ), typeof(ϕ), eltype(eltype(dirs)))
    return FeatureGaussian(SVector{length(μ),t}(μ), t(σ), t(ϕ), SVector{length(μ),t}[SVector{length(μ),t}(dir/norm(dir)) for dir in dirs], nodes, nodeidxs, type)
end

convert(::Type{FeatureGaussian{N,T,K}}, g::FeatureGaussian) where {N,T,K} = FeatureGaussian(SVector{N,T}(g.μ), T(g.σ), T(g.ϕ), Vector{SVector{N,T}}(g.dirs), g.nodes, g.nodeidxs, g.type)
convert(::Type{IsotropicGaussian{N,T}}, g::FeatureGaussian) where {N,T} = IsotropicGaussian(SVector{N,T}(g.μ), T(g.σ), T(g.ϕ), Vector{SVector{N,T}}(g.dirs))

function transform!(fg::FeatureGaussian, tform::AffineMap)
    fg.μ = tform(fg.μ)
    for (i,dir) in enumerate(fg.dirs)
        fg.dirs[i] = tform.linear * dir
    end
    return fg
end

# function (tform::AffineMap)(fg::FeatureGaussian)
#     new_fg = copy(fg)
#     transform!(new_fg, tform)
#     return new_fg
# end

"""
    feat = atoms_to_feature(atoms, nodeidxs, type, dirs=nothing)

Combines Gaussian distributions specified by `atoms` to create a feature repressented by a 
single Gaussian. The feature is centered at the center of mass of the distributions, and the
width `σ` of the feature is the average of the distance of individual features from the center 
added to their widths.

Geometric constraints for the feature can be specified by `dirs`. The returned feature has no
direction by default.
"""
function atoms_to_feature(gaussians::AbstractVector{<:AbstractIsotropicGaussian}, atoms, nodeidxs, type, dirs=nothing)
    if length(gaussians)==1
        μ = gaussians[1].μ
        σ = gaussians[1].σ
        ϕ = gaussians[1].ϕ
    else
        μ = center_of_mass(gaussians)
        ϕ = sum([a.ϕ for a in gaussians])
        σdiag = (sum([(a.σ^2 .+ a.μ.^2).*a.ϕ for a in gaussians]) ./ ϕ) .- μ.^2      # https://stats.stackexchange.com/questions/43159/how-to-calculate-pooled-variance-of-two-or-more-groups-given-known-group-varianc
        σ = prod(σdiag) ^ (1/length(σdiag))                                     # create spherical gaussian with the same volume as the anisotropic one represented by σdiag
    end
    if isnothing(dirs)
        dirs = typeof(μ)[]
    else
        dirs = [SVector{length(dir)}(dir) for dir in dirs]
    end

    return FeatureGaussian(μ, σ, ϕ, dirs, atoms, nodeidxs, type)
end
atoms_to_feature(mol::Union{UndirectedGraph,SubgraphView}, nodeidxs, σfun, ϕfun, args...
    ) = atoms_to_feature(MolGMM(mol, nodeidxs; σfun=σfun, ϕfun=ϕfun), [nodeattr(mol,i) for i in nodeidxs], nodeidxs, args...)
atoms_to_feature(gmm::AbstractSingleGMM, args...) = atoms_to_feature(gmm.gaussians, args...)
