"""
    tform = gogma_tform(X)

Converts the angle-axis rotation and translation parameters from GOGMA to an `AffineMap`.
"""
function gogma_tform(X)
    return AffineMap(GOGMA.rotmat(X[1],X[2],X[3]), SVector(X[4],X[5],X[6]))
end

"""
    score, tform, nevals = gogma_align_gmms(gmmfixed, gmmmoving, nsplits=2; atol=0.1, rtol=0.01, maxblocks=5e8, maxstagnant=Inf, maxevals=Inf)

Finds the optimal alignment between the two supplied GMMs a branch and bound method,
based on the [GOGMA alignment algorithm.](https://arxiv.org/abs/1603.00150)
"""
function gogma_align_gmms(gmm1::Union{MolGMM,FeatureMolGMM}, gmm2::Union{MolGMM,FeatureMolGMM}, nsplits=2; atol=0.1, rtol=0.01, maxblocks=5e8, maxstagnant=Inf, maxevals=Inf)
    res = GOGMA.gogma_align(gmm2.model, gmm1.model, nsplits, atol=atol, rtol=rtol, maxblocks=maxblocks, maxevals=maxevals, maxstagnant=maxstagnant, threads=threads)
    return res[1], gogma_tform(res[3]), res[4]
end

"""
    score, tform, nevals = rot_gogma_align_gmms(gmmfixed, gmmmoving, nsplits=2; atol=0.1, rtol=0.01, maxblocks=5e8, maxstagnant=Inf, maxevals=Inf)

Finds the optimal rotation between the two supplied GMMs a branch and bound method, for a given translation,
based on the [GOGMA alignment algorithm.](https://arxiv.org/abs/1603.00150)
"""
function rot_gogma_align_gmms(gmm1::Union{MolGMM,FeatureMolGMM}, gmm2::Union{MolGMM,FeatureMolGMM}, nsplits=2, trl=nothing; atol=0.1, rtol=0.01, maxblocks=5e8, maxstagnant=Inf, maxevals=Inf)
    res = GOGMA.rot_gogma_align(gmm2.model, gmm1.model, nsplits, trl, atol=atol, rtol=rtol, maxblocks=maxblocks, maxevals=maxevals, maxstagnant=maxstagnant)
    rx, ry, rz = res[3]
    R = GOGMA.rotmat(rx, ry, rz)
    if isnothing(trl)
        tform = AffineMap(R, @SVector zeros(size(gmm1,2)))
    else
        tform = AffineMap(R, trl)
    end
    return res[1], tform, res[4]
end

"""
    score, tform, nevals = tiv_gogma_align_gmms(gmmfixed, gmmmovingn, splits=2; atol=0.1, rtol=0.1, maxblocks=5e8, maxstagnant=Inf, maxevals=Inf)

Finds the optimal alignment between the two supplied GMMs using GOGMA and translation invariant vectors (TIVs),
inspired by [TIV-3DIV.](https://arxiv.org/abs/1812.11307). `cfixed` and `cmoving` define the ratio between numbers of TIVs used to represent
`gmmfixed` and `gmmmoving` respectively and the number of Gaussians they contain, during rotational alignment. 
"""
function tiv_gogma_align_gmms(gmmfixed::Union{MolGMM,FeatureMolGMM}, gmmmoving::Union{MolGMM,FeatureMolGMM}, cfixed=Inf, cmoving=Inf, nsplits=2; atol=0.1, rtol=0.01, maxblocks=5e8, maxstagnant=Inf, maxevals=Inf)
    res = GOGMA.tiv_gogma_align(gmmmoving.model, gmmfixed.model, cmoving, cfixed, nsplits, atol=atol, rtol=rtol, maxblocks=maxblocks, maxevals=maxevals, maxstagnant=maxstagnant)
    return res[1], gogma_tform(res[3]), res[4]
end
