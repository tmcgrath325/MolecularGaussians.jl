import Base: eltype, convert
import GaussianMixtureAlignment: combine, dims, numbertype

mutable struct AtomGaussian{N,T<:Real} <: AbstractIsotropicGaussian{N,T}
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
convert(::Type{IsotropicGaussian{N,T}}, g::AtomGaussian) where {N,T} = IsotropicGaussian(SVector{N,T}(g.μ), T(g.σ), T(g.ϕ), Vector{SVector{N,T}}(g.dirs))

function transform!(ag::AtomGaussian, tform::AffineMap)
    ag.μ = tform(ag.μ)
    for (i,dir) in enumerate(ag.dirs)
        ag.dirs[i] = tform.linear * dir
    end
    return ag
end

# function (tform::AffineMap)(ag::AtomGaussian)
#     new_ag = copy(ag)
#     transform!(new_ag, tform)
#     return new_ag
# end

