"""
    com = center_of_mass(gmm)

Returns the center of mass of `gmm`, where its first order moments are equal to 0.
"""
function center_of_mass(positions::AbstractMatrix{<:Real}, weights=ones(eltype(positions),size(positions,2)))
    return [sum(hcat([weights[i]*positions[:,i] for i=1:size(positions,2)]...), dims=2)...]/sum(weights)
end
center_of_mass(gaussians::AbstractVector{IsotropicGaussian{T,N}} where T where N, tforms=fill(identity, length(gaussians)), weights=[g.ϕ for g in gaussians]) = 
    return center_of_mass(hcat([tforms[i](g.μ) for (i,g) in enumerate(gaussians)]...), weights)
center_of_mass(gmm::IsotropicGMM, tform::Union{typeof(identity),AffineMap}=identity) = center_of_mass(gmm.gaussians, fill(tform, length(gmm)))
center_of_mass(gmm::MolGMM, tform::Union{typeof(identity),AffineMap}=identity) = center_of_mass(gmm.model, tform)

""" 
    m = second_moment(gmm, center, dim1, dim2)

Returns the second order moment of `gmm`
"""
function mass_matrix(positions::AbstractMatrix{<:Real}, 
                     widths=zeros(eltype(positions),size(positions,2)), 
                     weights=ones(eltype(positions),size(positions,2)), 
                     center=center_of_mass(positions,weights))
    t = eltype(positions)
    npts = size(positions,2)
    M = fill(zero(t), 3, 3)
    dists = [positions[:,i].-center for i=1:npts]
    # diagonal terms
    for i=1:3
        for j=1:npts
            M[i,i] += weights[j] * (dists[j][i]^2 + widths[j]^2)
        end
    end
    # off-diagonal terms
    for i=1:3
        for j=i+1:3
            for k=1:npts
                inc = weights[k] * dists[k][i] * dists[k][j]
                M[i,j] += inc
                M[j,i] += inc
            end
        end
    end
    return M
end
mass_matrix(gmm::IsotropicGMM, tform::Union{typeof(identity),AffineMap}=identity, center=center_of_mass(gmm, tform)) =
    mass_matrix(hcat([tform(g.μ) for g in gmm.gaussians]...), [g.σ for g in gmm.gaussians], [g.ϕ for g in gmm.gaussians], center)

mass_matrix(gmm::MolGMM, tform=identity, com=center_of_mass(gmm,tform)) = mass_matrix(gmm.model, tform, com)

"""
    tforms = inertial_transforms(gmm)

Returns 4 transformations to put `gmm` in an inertial frame. That is, the mass matrix of the GMM
is made diagonal, and the GMM center of mass is made the origin.
"""
function inertial_transforms(positions::AbstractMatrix{<:Real},
                             widths=zeros(eltype(positions),size(positions,2)), 
                             weights=ones(eltype(positions),size(positions,2)); 
                             invert = false)
    com = center_of_mass(positions, weights)
    massmat = mass_matrix(positions, widths, weights, com)
    evals, evecs = eigvals(massmat), eigvecs(massmat)
    eorder = sortperm(evals, rev=true)
    R1 = zeros(eltype(massmat),size(massmat))
    for i=1:size(massmat,1)
        R1[:,i] = evecs[:,eorder[i]]
    end
    R1 = SMatrix{3,3}(R1)
    R2 = R1 * @SMatrix [1 0 0; 0 -1 0; 0 0 -1]
    R3 = R1 * @SMatrix [-1 0 0; 0 1 0; 0 0 -1]
    R4 = R1 * @SMatrix [-1 0 0; 0 -1 0; 0 0 1]
    if invert
        return AffineMap(R1, com), AffineMap(R2, com), AffineMap(R3, com), AffineMap(R4, com)
    else
        return AffineMap(transpose(R1), -com), AffineMap(transpose(R2), -com), AffineMap(transpose(R3), -com), AffineMap(transpose(R4), -com)
    end
end
inertial_transforms(gmm::IsotropicGMM; invert=false) = 
    return inertial_transforms(hcat([g.μ for g in gmm.gaussians]...), [g.σ for g in gmm.gaussians], [g.ϕ for g in gmm.gaussians]; invert=invert)
inertial_transforms(gmm::MolGMM; invert=false) = inertial_transforms(gmm.model; invert=invert)

"""
    score, tform, nevals = rocs_align_gmms(gmmfixed, gmmmoving; maxevals=1000)

Finds the optimal alignment between the two supplied GMMs using steric multipoles,
based on the [ROCS alignment algorithm.](https://docs.eyesopen.com/applications/rocs/theory/shape_shape.html)
"""
function rocs_align_gmms(gmmfixed::MolGMM, gmmmoving::MolGMM; maxevals=1000)
    # align both GMMs to their inertial frames
    tformfixed = inertial_transforms(gmmfixed)[1]  # Only need one alignment for the fixed GMM
    tformsmoving = inertial_transforms(gmmmoving)    # Take all 4 rotations for the GMM to be aligned

    # perform local alignment starting at each inertial rotation for the "moving" GMM
    results = []
    for i=1:length(tformsmoving)
        push!(results, local_align_gmms(gmmfixed, gmmmoving, tformfixed=tformfixed, tformmoving=tformsmoving[i], maxevals=maxevals))
    end 

    # Apply the inverse of `tformfixed` to the optimized transformation 
    minoverlap, mindex = findmin([r[1] for r in results])
    tformmoving = results[mindex][2]
    nevals = sum([r[3] for r in results])
    R = transpose(tformfixed.linear) * tformmoving.linear
    tr = -tformfixed.translation + tformmoving.translation
    return minoverlap, AffineMap(R, tr), nevals
end

rocs_align_gmms(molfixed::Union{UndirectedGraph, SubgraphView}, molmoving::Union{UndirectedGraph, SubgraphView}, σfun=ones, ϕfun=ones; maxevals=1000) = 
    rocs_align_gmms(MolGMM(molfixed, σfun, ϕfun), MolGMM(molmoving, σfun, ϕfun), maxevals=maxevals)