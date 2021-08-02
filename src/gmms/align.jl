"""
    objval = alignment_objective(X, gmmfixed, gmmmoving)
    X = (rx, ry, rz, tx, ty, tz)

Calculates the "overlap" between two three-dimensional Gaussian Mixture Models, `gmmfixed` 
and `gmmmoving`, with a rigid transformation defined by vector X applied to `gmmmoving`.

Within the vector `X`, `tx`, `ty`, and `tz` correspond to the components of a translations in 
three-dimensional space, and `rx`, `ry`, and `rz` are the vector components of a 
[angle-axis representation](https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation)
""" 
function alignment_objective(X, gmmfixed::Union{MolGMM,FeatureMolGMM}, gmmmoving::Union{MolGMM,FeatureMolGMM})
    t = promote_type(eltype(gmmfixed), eltype(gmmmoving))
    if sum(abs2, X[1:3]) > t(π)
        return typemax(t)
    end
    tform = gogma_tform(X)
    return -gmm_overlap(gmmfixed, gmmmoving, tform2=tform)
end
