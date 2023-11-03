

"""
    pharmgmm = commonfeature_pharmacophore(mols, active, min_features, max_features)

Creates a pharmacophore model from a set of `mols` which consists of the most commonly shared features
following alignment. The `active` vector indicates which molecules are active and which are inactive.
"""
function commonfeature_pharmacophore(template::PharmacophoreGMM{N,T,K}, confs::AbstractVector{<:AbstractVector{<:PharmacophoreGMM{N,T,K}}}, active::AbstractVector{Bool}; 
                                     max_totalfeatures::Int=3, 
                                     min_totalfeatures::Int=1, 
                                     max_excluded_volumes::Int=2,
                                     min_features::Union{Dict{Symbol,Int},Nothing}=nothing, 
                                     max_features::Union{Dict{Symbol,Int},Nothing}=nothing) where {N,T,K}
    @assert min_totalfeatures <= max_totalfeatures
    min_features = isnothing(min_features) ? Dict{Symbol,Int}() : min_features
    max_features = isnothing(max_features) ? Dict{Symbol,Int}() : max_features
    @assert sum(values(min_features)) <= max_totalfeatures
    # align all mols to the template
    mols_with_tform = [align_conformers(cs, template) for cs in confs]
    aligned_mols = [tformgraph(mtf[2](mtf[1]), mtf[2]) for mtf in mols_with_tform]
    # group features by type
    grouped_features = Dict{K, Tuple{Vector{IsotropicGaussian{N,T}},Vector{Bool}}}()
    for (i,mol) in enumerate(aligned_mols)
        for key in keys(mol.featuremaps)
            if !haskey(grouped_features, key)
                grouped_features[key] = (Vector{IsotropicGaussian{N,T}}(), Vector{Bool}())
            end
            gmm = mol.gmms[key]
            push!(grouped_features[key][1], gmm.gaussians...)
            push!(grouped_features[key][2], fill(active[i], length(gmm.gaussians))...)
        end
    end
    # Calculate and rank overlap scores for each feature in the template
    template_olaps = Dict{K, Vector{Float64}}()
    template_rankings = Dict{K, Vector{Int}}()
    for (key, (features, active)) in grouped_features
        olaps = template_overlap_scores(template.gmms[key].gaussians, features, active)
        push!(template_olaps, key => olaps)
        push!(template_rankings, key => sortperm(template_olaps[key]; rev=true))
    end
    # initialize model
    pharmacophore = IsotropicMultiGMM{N,T,K}(Dict{K, IsotropicGMM{N,T}}())
    model_gmms = pharmacophore.gmms
    n_features_added = 0
    feature_keys = collect(keys(grouped_features))
    for key in feature_keys
        model_gmms[key] = IsotropicGMM{N,T}([])
    end
    # satisfy minimum feature requirements
    for key in keys(min_features)
        rankings = template_rankings[key]
        for i=1:min_features[key]
            push!(model_gmms[key].gaussians, template.gmms[key].gaussians[rankings[i]])
            n_features_added += 1
        end
    end
    # add features until maximum feature constraints are met
    template_feat_counts = [length(template.gmms[key].gaussians) for key in feature_keys]
    while n_features_added < max_totalfeatures
        next_idxs = [length(model_gmms[key].gaussians)+1 for key in feature_keys]
        best_scores = [idx > temp_max ? -Inf : template_olaps[key][template_rankings[key][idx]] for (key,temp_max,idx) in zip(feature_keys, template_feat_counts, next_idxs)]
        next_key_idx = argmax(best_scores)
        next_key = feature_keys[next_key_idx]
        while haskey(max_features, next_key) && (max_features[next_key] < next_idxs[next_key_idx]) && (maximum(best_scores) != -Inf)
            best_scores[next_key_idx] = -Inf
            next_key_idx = argmax(best_scores)
            next_key = feature_keys[next_key_idx]
        end
        best_scores[next_key_idx] == -Inf && break
        feature_to_add = template.gmms[next_key][next_idxs[next_key_idx]]
        push!(model_gmms[feature_keys[next_key_idx]].gaussians, feature_to_add)
        n_features_added += 1
    end
    # add excluded volume features
    if (max_excluded_volumes > 0) && (n_features_added < max_totalfeatures)
        push!(model_gmms, :Volume => IsotropicGMM{N,T}([]))
        volume_features = reduce(vcat, [feats[1] for feats in values(grouped_features)])
        volume_active = reduce(vcat, [feats[2] for feats in values(grouped_features)])
        dmat = feature_distance_matrix(volume_features)
        omat = feature_overlap_matrix(volume_features, dmat, volume_active)
        while (length(model_gmms[:Volume]) < max_excluded_volumes) && (n_features_added < max_totalfeatures)
            feat_overlap_scores = sum(omat, dims=2)
            idx = argmin(feat_overlap_scores)
            push!(model_gmms[:Volume].gaussians, volume_features[idx])
            n_features_added += 1
            omat[idx,:] .= 0
            omat[:,idx] .= 0
        end
    end
    # remove empty sets of features
    for key in feature_keys
        if length(pharmacophore.gmms[key].gaussians) == 0
            delete!(pharmacophore.gmms, key)
        end
    end
    return pharmacophore
end

function template_overlap_scores(template_features::AbstractVector{<:G}, features::AbstractVector{<:G}, active) where {N,T,G<:IsotropicGaussian{N,T}}
    scores = fill(zero(T), length(template_features))
    for (i,t) in enumerate(template_features)
        for (j,f) in enumerate(features)
            scores[i] += active[j] ? (overlap(f,t) / sum(active)) : (-overlap(f,t) / (length(active)-sum(active)))
        end
    end
    return scores
end

function distances_from_template_features(template_features::AbstractVector{<:G}, features::AbstractVector{<:G}) where {N,T,G<:IsotropicGaussian{N,T}}
    tlen = length(template_features)
    flen = length(features)
    dmat = zeros(tlen, flen)
    for i=1:tlen
        for j=1:flen
            dmat[i,j] = norm(template_features[i].μ .- features[j].μ)
            dmat[j,i] = dmat[i,j]
        end
    end
    return dmat
end

function feature_distance_matrix(features)
    len = length(features)
    dmat = zeros(len, len)
    for i=1:len
        for j=i+1:len
            dmat[i,j] = norm(features[i].μ .- features[j].μ)
            dmat[j,i] = dmat[i,j]
        end
    end
    return dmat
end

function feature_overlap_matrix(features, dmat, active)
    len = length(features)
    omat = zeros(len, len)
    for i=1:len
        for j=i+1:len
            omat[i,j] = (active[i] && active[j]) ?  (overlap(dmat[i,j], features[i].σ, features[j].σ, features[i].ϕ, features[j].ϕ) / sum(active)) : 
                                                    (overlap(dmat[i,j], features[i].σ, features[j].σ, features[i].ϕ, features[j].ϕ) / (length(active)-sum(active)))
            omat[j,i] = omat[i,j]
        end
    end
    return omat
end

