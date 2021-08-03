"""
    intg = gm_product_integral(gm1, gm2; tform1, tform2) 

Computes the integral of the product of two spherical Gaussian distributions, `gm1` and `gm2`, across 
n-dimensional space.

If provided, rigid transformations `tform1` and `tform2` will be applied to `gm1` and `gm2` respectively
"""
function gm_product_integral(gm1::IsotropicGaussian, gm2::IsotropicGaussian, tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity)
    if length(gm1.dirs) * length(gm2.dirs) == 0
        dirdot = one(promote_type(eltype(gm1),eltype(gm2)))
    else
        dirdot = -one(promote_type(eltype(gm1),eltype(gm2)))
        for dir1 in gm1.dirs
            for dir2 in gm2.dirs
                if typeof(tform1) == typeof(identity)
                    tdir1 = dir1
                else
                    tdir1 = tform1.linear*dir1
                end
                if typeof(tform2) == typeof(identity)
                    tdir2 = dir2
                else
                    tdir2 = tform2.linear*dir2
                end
                dirdot = max(dirdot, dot(tdir1,tdir2))
            end
        end
    end
    return -GaussianMixtureAlignment.objectivefun(sum(abs2, tform1(gm1.μ)-tform2(gm2.μ)), 
                               gm1.σ^2 + gm2.σ^2, 
                               gm1.ϕ * gm2.ϕ,
                               dirdot)
end

"""
    ovrlp = gmm_overlap(gmm1, gmm2; tform1=identity, tform2=identity)

    Calculates the "overlap" between two Gaussian Mixture Models, which is the sum of
    the integrals of products of each pairwise combination of distributions within the
    two models.
"""
function gmm_overlap(gmm1::IsotropicGMM, gmm2::IsotropicGMM; tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity)
    val = zero(promote_type(eltype(gmm1), eltype(gmm2)))
    for gm1 in gmm1.gaussians
        for gm2 in gmm2.gaussians
            val = val + gm_product_integral(gm1, gm2, tform1, tform2)
        end
    end
    return val
end

function gmm_overlap(mgmm1::MultiGMM, mgmm2::MultiGMM, ks=keys(mgmm1.gmms) ∪ keys(mgmm2.gmms); tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity)
    o = zero(promote_type(eltype(mgmm1),eltype(mgmm2)))
    for key in ks
        if (haskey(mgmm1.gmms, key) && haskey(mgmm2.gmms, key))
            o += gmm_overlap(mgmm1.gmms[key], mgmm2.gmms[key], tform1=tform1, tform2=tform2)
        end
    end
    return o
end

gmm_overlap(gmm1::Union{MolGMM,FeatureMolGMM}, gmm2::Union{MolGMM,FeatureMolGMM}; tform1=identity, tform2=identity) = gmm_overlap(gmm1.model, gmm2.model, tform1=tform1, tform2=tform2)

"""
    l2dist = gmm_distance(model1, model2)

Calculates the L2 distance between two Gaussian Models made up of spherical Gaussian distributions.
"""
function gmm_distance(gmm1::IsotropicGMM, gmm2::IsotropicGMM; tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity)
    return gmm_overlap(gmm1, gmm1, tform1=identity, tform2=identity) + 
           + gmm_overlap(gmm2, gmm2, tform1=identity, tform2=identity) +
           -2*gmm_overlap(gmm1, gmm2, tform1=tform1, tform2=tform2)
end

function gmm_distance(mgmm1::MultiGMM, mgmm2::MultiGMM, ks=keys(mgmm1.gmms) ∪ keys(mgmm2.gmms); tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity)
    return gmm_overlap(mgmm1, mgmm1, ks; tform1=identity, tform2=identity) + 
           + gmm_overlap(mgmm2, mgmm2, ks; tform1=identity, tform2=identity) +
           -2*gmm_overlap(mgmm1, mgmm2, ks; tform1=tform1, tform2=tform2)
end

gmm_distance(molgmm1::Union{MolGMM,FeatureMolGMM}, molgmm2::Union{MolGMM,FeatureMolGMM}; tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity) =
    gmm_distance(molgmm1.model, molgmm2.model; tform1=tform1, tform2=tform2)

gmm_distance(fgmm1::FeatureMolGMM, fgmm2::FeatureMolGMM, features::AbstractVector{<:Symbol}; tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity) =
    gmm_distance(fgmm1.model, fgmm2.model, features; tform1=tform1, tform2=tform2)

gmm_distance(mol1::Union{UndirectedGraph,SubgraphView}, mol2::Union{UndirectedGraph,SubgraphView}, σfun=ones, ϕfun=ones; tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity) = 
    gmm_distance(MolGMM(mol1, σfun, ϕfun), MolGMM(mol2, σfun, ϕfun), tform1=tform1, tform2=tform2)


function gmm_tanimoto(molgmm1::MolGMM, molgmm2::MolGMM; tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity)
    o = gmm_overlap(molgmm1, molgmm2, tform1=tform1, tform2=tform2)
    return o / (gmm_overlap(molgmm1, molgmm1, tform1=identity, tform2=identity) + gmm_overlap(molgmm2, molgmm2, tform1=identity, tform2=identity) - o)
end

function gmm_tanimoto(fgmm1::FeatureMolGMM, fgmm2::FeatureMolGMM; tform1::Union{typeof(identity),AffineMap}=identity, tform2::Union{typeof(identity),AffineMap}=identity)
    o = zero(promote_type(eltype(fgmm1),eltype(fgmm2)))
    self = zero(promote_type(eltype(fgmm1),eltype(fgmm2)))
    for key in keys(fgmm1.gmms) ∪ keys(fgmm2.gmms)
        if (haskey(fgmm1.gmms, key) && haskey(fgmm2.gmms, key))
            o += gmm_overlap(fgmm1.gmms[key], fgmm2.gmms[key], tform1=tform1, tform2=tform2)
        end
        if haskey(fgmm1.gmms, key)
            self += gmm_overlap(fgmm1.gmms[key], fgmm1.gmms[key], tform1=identity, tform2=identity)
        end
        if haskey(fgmm2.gmms, key)
            self += gmm_overlap(fgmm2.gmms[key], fgmm2.gmms[key], tform1=identity, tform2=identity)
        end
    end
    return o / (self-o)
end
