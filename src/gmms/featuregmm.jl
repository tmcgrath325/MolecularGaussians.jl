
mutable struct FeatureGMM{N,T<:Real,K} <: AbstractIsotropicGMM{N,T}
    gaussians::Vector{FeatureGaussian{N,T,K}}
    nodes::Vector{Set{Int}}
    type::K
end

eltype(::Type{FeatureGMM{N,T,K}}) where {N,T,K} = FeatureGaussian{N,T,K}

convert(::Type{FeatureGMM{N,T,K}}, f::FeatureGMM) where {N,T,K} = FeatureGMM([convert(FeatureGaussian{N,T,K},g) for g in f.gaussians], f.nodes, f.type)
convert(::Type{IsotropicGMM{N,T}}, f::FeatureGMM) where {N,T} = IsotropicGMM([convert(IsotropicGaussian{N,T},g) for g in f.gaussians])

function combine(featgmm1::FeatureGMM, featgmm2::FeatureGMM)
    if dims(featgmm1) != dims(featgmm2)
        throw(ArgumentError("GMMs must have the same dimensionality"))
    end
    if featgmm1.type != featgmm2.type
        throw(ArgumentError("Each FeatureGMM must have the same feature type"))
    end
   
    t = promote_type(eltype(featgmm2), eltype(featgmm2))
    gaussians = t[vcat(featgmm1.gaussians, featgmm2.gaussians)...]

    nodes = append(featgmm1.nodes, featgmm2.nodes)

    return FeatureGMM(gaussians, nodes, featgmm1.type)
end


# rigid transformation

function transform!(featgmm::FeatureGMM, tform::AffineMap)
    for fg in featgmm.gaussians
        transform!(fg, tform)
    end
    return featgmm
end

# function (tform::AffineMap)(featgmm::FeatureGMM)
#     new_featgmm = deepcopy(featgmm)
#     transform!(new_featgmm, tform)
#     return new_featgmm
# end


# bond rotations

function bondrotate!(featgmm::FeatureGMM, rotsubgraph::RotableSubgraph, tform)
    for (j,fg) in enumerate(featgmm.gaussians)
        if prod([i ∈ rotsubgraph.nodeidxs for i in fg.nodeidxs])
            featgmm.gaussians[j] = tform(fg)
        end
    end
    return featgmm
end

function bondrotate(featgmm::FeatureGMM{N,T,K}, rotsubgraph::RotableSubgraph, tform) where {N,T,K}
    n = N
    t = promote_type(T, eltype(tform.linear), eltype(tform.translation))
    k = K

    gaussians = FeatureGaussian{n,t,k}[]
    for (j,fg) in enumerate(featgmm.gaussians)
        if prod([i ∈ rotsubgraph.nodeidxs for i in fg.nodeidxs])
            push!(gaussians, tform(fg))
        else
            push!(gaussians, fg)
        end
    end
    return FeatureGMM(gaussians, featgmm.nodes, featgmm.type)
end

# function bondsrotate!(featgmm::FeatureGMM, rotbondidxs, angles)
#     for (i,ridx) in enumerate(rotbondidxs)
#         bondrotate!(featgmm, ridx, angles[i])
#     end
#     return featgmm
# end

# function bondsrotate(featgmm::FeatureGMM, rotbondidxs, angles)
#     for (i,ridx) in enumerate(rotbondidxs)
#         featgmm = bondrotate(featgmm, ridx, angles[i])
#     end
#     return featgmm
# end


# descriptive display

Base.show(io::IO, featgmm::FeatureGMM) = println(io,
    "FeatureGMM from molecule with formula $(typeof(featgmm.graph)<:GraphMol ? molecularformula(featgmm.graph) : molecularformula(featgmm.graph.graph))",
    " with $(length(featgmm)) Gaussians."
)