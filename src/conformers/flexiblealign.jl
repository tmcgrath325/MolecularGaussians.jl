## Bridge a PharmacophoreGMM to GaussianMixtureAlignment's articulated-alignment interface, so
## a molecule's rotatable bonds act as the joints of a flexible branch-and-bound alignment. The
## bond geometry PharmacophoreGMM already stores (`axes`, `origins`, `bondtogaussians`,
## `bondtobonds`) supplies every quantity the interface needs; `bondrotate` is the forward
## kinematics. Rotatable bonds are ordered so each precedes the bonds downstream of it, which is
## the ordering the interface requires.

njoints(pgmm::PharmacophoreGMM) = length(pgmm.axes)
joint_axis(pgmm::PharmacophoreGMM, b) = normalize(pgmm.axes[b])
joint_origin(pgmm::PharmacophoreGMM, b) = pgmm.origins[b]
joint_features(pgmm::PharmacophoreGMM, b) = pgmm.bondtogaussians[b]
joint_children(pgmm::PharmacophoreGMM, b) = pgmm.bondtobonds[b]
flex(pgmm::PharmacophoreGMM, φ) = bondrotate(pgmm, φ, eachindex(pgmm.axes))

"""
    flexible_align(x::PharmacophoreGMM, y; kwargs...)

Globally align molecule model `x` onto `y`, optimizing the rigid rotation and translation
together with a rotation about each of `x`'s rotatable bonds. Returns a
`FlexibleAlignmentResult`; `aligned(result)` is the posed, flexed `x`.

This is a thin wrapper over `flex_gogma_align`; keyword arguments (tolerances, iteration limits,
`interactions`) are forwarded to it. Unlike enumerating discrete conformers and aligning each
rigidly, the bond angles are searched continuously alongside the rigid transform.
"""
flexible_align(x::PharmacophoreGMM, y::AbstractGMM; kwargs...) = flex_gogma_align(x, y; kwargs...)
