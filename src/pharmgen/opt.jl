function feature_points(template::AbstractIsotropicMultiGMM{N,T}) where {N,T}
    positions = zeros(T,sum(length, values(template.gmms)),N)
    idx = 1
    for (k, gmm) in template.gmms
        for g in gmm.gaussians
            positions[idx,:] = g.μ
            idx += 1
        end
    end
    return positions
end

# input assumes that the template has been aligned to the coordinate axes
function grid_points(template::AbstractIsotropicMultiGMM{N,T}, sizecoef=1.5, spacing=1.0) where {N,T}
    feat_positions = feature_points(template)
    max_pos = maximum(feat_positions, dims=1)
    min_pos = minimum(feat_positions, dims=1)
    center = vec(max_pos .+ min_pos) ./ 2
    gridsize = (max_pos .- min_pos) .* sizecoef
    return grid_points(center, gridsize, spacing)
end

function grid_points(center::AbstractVector{T}, gridsize, spacing=1.0) where T
    N = length(center)
    @assert length(gridsize) == N
    excess = gridsize .% spacing
    if excess != 0
        gridsize = gridsize .- excess
    end
    npoints = tuple(map(Int, (gridsize ./ spacing) .+ 1)...)
    factors = [1, cumprod(npoints[1:end-1])...]
    positions = zeros(T, prod(npoints), N)
    for I in CartesianIndices(npoints)
        idx = Tuple(I)
        i = sum(factors .* (idx .-1)) + 1
        positions[i,:] = center .+ (idx .- 1) .* spacing
    end
    return positions
end

function labeled_points(positions::AbstractMatrix, labels::AbstractVector)
    n = length(labels)
    m = size(positions)[1]
    lpositions = repeat(positions, outer=n)
    plabels = repeat(labels; inner=m)
    return lpositions, plabels
end

function activity_prediction(molgmm::G, modelgmm::L, conc::Real, max_activity::Real) where {G,L}
    overlap_score = overlap(molgmm, modelgmm)
    return max_activity * conc / (overlap_score + conc)
end

function activity_prediction(molgmm::G, modelgmm::L, concs::AbstractVector{<:Real}, max_activity::Real) where {G,L}
    return [activity_prediction(molgmm, modelgmm, conc, max_activity) for conc in concs]
end

function activity_prediction(models::AbstractVector{<:G}, concs::AbstractVector{<:AbstractVector{<:Real}}, modelgmm::L, max_activity::Real) where {G,L}
    return [activity_prediction(mol, modelgmm, conc, max_activity) for (mol, conc) in zip(models, concs)]
end

function prediction_sq_error(models::AbstractVector{<:G}, concs::AbstractVector{<:AbstractVector{<:Real}}, modelgmm::L, activities::AbstractVector{<:AbstractVector{<:Real}}) where {G,L}
    pred = activity_prediction(models, concs, modelgmm, maximum(maximum, activities))
    return sum((collect(Iterators.flatten(activities)) .- collect(Iterators.flatten(pred))).^2)
end

function lasso_penalty(weights::AbstractVector{<:Real}, λ::Real)
    return λ * sum(abs.(weights))
end

function features_from_labeled_points(positions::AbstractMatrix{T}, labels::AbstractVector{<:L}, weights::AbstractVector{<:S}=ones(length(labels)), widths=ones(length(labels))) where {T,L,S}
    N = size(positions)[2]
    @assert size(positions)[1] > 0 && N > 0
    gmms = Dict{eltype(labels), IsotropicGMM{N,S}}()
    for (i,(label, σ, ϕ)) in enumerate(zip(labels, widths, weights))
        μ = positions[i,:]
        if !haskey(gmms, label)
            gmms[label] = IsotropicGMM{N,S}([])
        end
        if ϕ !== 0 # ignore zero-weighted features
            push!(gmms[label].gaussians, IsotropicGaussian{N,S}(μ, σ, ϕ))
        end
    end
    return IsotropicMultiGMM(gmms)
end

function weights_cost_fun(modelgmm::G, weights, models::AbstractVector{<:G}, concs::AbstractVector{<:AbstractVector{<:Real}}, activities::AbstractVector{<:AbstractVector{<:Real}}, λ::Real) where G
    return prediction_sq_error(models, concs, modelgmm, activities) + lasso_penalty(weights, λ)
end
weights_cost_fun(positions, labels, weights, models, concs, activities, λ=0.0) = weights_cost_fun(features_from_labeled_points(positions, labels, weights), weights, models, concs, activities, λ)

function relaxation_term(weights, σ)
    return sum([(w + σ)^2 for w in weights])
end

function relaxed_weights_cost_fun(modelgmm::G, weights, models::AbstractVector{<:G}, concs::AbstractVector{<:AbstractVector{<:Real}}, activities::AbstractVector{<:AbstractVector{<:Real}}, λ::Real, σ::Real) where G
    return prediction_sq_error(models, concs, modelgmm, activities) + lasso_penalty(weights, λ) + relaxation_term(weights, σ)   
end
relaxed_weights_cost_fun(positions, labels, weights, models, concs, activities, λ, σ) = relaxed_weights_cost_fun(features_from_labeled_points(positions, labels, weights), weights, models, concs, activities, λ, σ)


function alignment_cost_fun(modelgmm::G, models::AbstractVector{<:G}) where G
    return sum([tanimoto(modelgmm, model) for model in models])
end
alignment_cost_fun(positions, labels, weights, models) = alignment_cost_fun(features_from_labeled_points(positions, labels, weights), models)

function optimize_weights_from_aligned_gmms(models::AbstractVector{<:G}, concs::AbstractVector{<:AbstractVector{<:Real}}, positions, labels, activities, λ) where {G<:AbstractGMM}
    @assert length(positions) === length(labels)
    init_weights = zeros(length(positions))
    cost_fun = (w) -> weights_cost_fun(positions, labels, w, models, concs, activities, λ)
    res = optimize(cost_fun, init_weights, LBFGS(), Optim.Options(iterations=100), autodiff=:forward)
    return res
end

function weights_relaxation(models::AbstractVector{<:G}, concs::AbstractVector{<:AbstractVector{<:Real}}, positions, labels, activities, λ, σ, α, its) where G
    @assert size(positions,1) === length(labels)
    init_weights = zeros(length(positions))
    w_i = init_weights
    σ_i = σ
    res = nothing
    print("Relaxation start")
    for i=1:3
        print("\r")
        print("Relaxation Iteration ", i, ":\t")
        cost_fun = (w) -> relaxed_weights_cost_fun(positions, labels, w, models, concs, activities, λ, σ_i)
        res = optimize(cost_fun, init_weights, LBFGS(), Optim.Options(iterations=100), autodiff=:forward)
        σ_i = σ_i * α
        w_i = res.minimizer
    end
    return res
end

function combined_cost_fun(tforms::AbstractVector{<:Real}, weights::AbstractVector{<:Real}, models::AbstractVector{<:G}, concs::AbstractVector{<:AbstractVector{<:Real}}, positions, labels, activities, β=1.0, λ=0.0) where {G<:AbstractGMM}
    am_tforms = [GaussianMixtureAlignment.AffineMap((tforms[(i-1)*6+1:i*6]...,)) for i=1:length(models)]
    tformed_models = [tform(model) for (tform, model) in zip(am_tforms, models)]
    modelgmm = features_from_labeled_points(positions, labels, weights)
    weights_cost = weights_cost_fun(modelgmm, weights, tformed_models, concs, activities, λ)
    alignment_cost = alignment_cost_fun(modelgmm, tformed_models)
    return weights_cost + β * alignment_cost
end

function optimize_combined_cost(conformers::AbstractVector{<:AbstractVector{<:G}}, concs::AbstractVector{<:AbstractVector{<:Real}}, positions, labels, activities, β=1.0, λ=0.0, σ=1.0, α=0.1; iterations = 100) where {G<:AbstractGMM}
    conf_idxs = [1 for mol in conformers]
    tforms = zeros(6 * length(conformers))
    models = [GaussianMixtureAlignment.AffineMap((tforms[(i-1)*6+1:i*6]...,))(molconfs[idx]) for (i, (idx, molconfs)) in enumerate(zip(conf_idxs, conformers))]
    weights = zeros(length(positions))
    
    for it=1:iterations
        print("Iteration ", it, ":\t")
        print("optimizing weights... ")
        # w_cost_fun = (w) -> weights_cost_fun(features_from_labeled_points(positions, labels, w), w, models, concs, activities, λ)
        # optimize weights
        # weights_res = optimize(w_cost_fun, weights, LBFGS(), Optim.Options(iterations=100), autodiff=:forward)
        # weights = weights_res.minimizer
        weights_res = weights_relaxation(models, concs, positions, labels, activities, λ, σ, α, 3)
        # update model from new weights
        modelgmm = features_from_labeled_points(positions, labels, weights)

        # optimize tforms and conformer choice
        for (i,molconfs) in enumerate(conformers)
            bestconf, tform, idx, min = align_conformers(molconfs, modelgmm; alignfun = local_align)
            conf_idxs[i] = idx
            for j=1:6
                tforms[(i-1)*6+j] = tform[j]
            end
            models[i] = bestconf
        end
        # combined optimization
        println("perfoming combined optimization... ")
        c_cost_fun = (p) -> combined_cost_fun(p[1:length(tforms)], p[length(tforms)+1:end], models, concs, positions, labels, activities, β, λ)
        res = optimize(c_cost_fun, [tforms; weights], LBFGS(), Optim.Options(iterations=100), autodiff=:forward)
        minimizer = res.minimizer
        tforms = minimizer[1:length(tforms)]
        weights = minimizer[length(tforms)+1:end]

        #update models with new tforms and weights
        modelgmm = features_from_labeled_points(positions, labels, weights)
        models = [GaussianMixtureAlignment.AffineMap((tforms[(i-1)*6+1:6*i]...,))(molconfs[idx]) for (i,(idx, molconfs)) in enumerate(zip(conf_idxs, conformers))]
        @show conf_idxs
        @show it, res.minimum, weights_res.minimum
    end
    return conf_idxs, tforms, weights, models, modelgmm
end