"""
MolecularGaussians.jl
===========================

MolecularGaussians.jl is a package used to align molecules and their pharmacophore features by modeling them as Gaussian mixture models.
It makes use of GaussianMixtureAlignment.jl to compute alignments, overlap, and distances between molecules.

REPL help
=========

? followed by an algorith or constructor name will print help to the terminal. See: \n
    \t?MolGMM \n
    \t?PharmacophoreGMM \n
    \t?gogma_align \n
    \t?tiv_gogma_align \n
    \t?rocs_align \n
"""
module MolecularGaussians

using LinearAlgebra

using StaticArrays
using CoordinateTransformations
using Rotations

using MolecularGraph
using Graphs

using GaussianMixtureAlignment
using GaussianMixtureAlignment: AbstractGaussian, AbstractSingleGMM, AbstractMultiGMM, AbstractGMM
using GaussianMixtureAlignment: AbstractIsotropicGaussian, AbstractIsotropicGMM, AbstractIsotropicMultiGMM
using GaussianMixtureAlignment: IsotropicGaussian, IsotropicGMM, IsotropicMultiGMM
using GaussianMixtureAlignment: centroid
using GaussianMixtureAlignment: local_align, rocs_align, gogma_align, rot_gogma_align, tiv_gogma_align, overlap, distance, tanimoto

export local_align, gogma_align, tiv_gogma_align, overlap, distance, tanimoto
export AtomGaussian, MolGMM, PharmacophoreGMM
export vdwradii, vdwradii!
export partialcharges, partialcharges!
export pharmfeatures, pharmfeatures!
export inertial_transforms, rocs_align
export affinetransform

export AtomType, FeatureDef

export moldisplay, pharmacophoregmmdisplay

using MakieCore
using GeometryBasics
using Colors

include("utils.jl")

include("features/featuredef.jl")
include("features/parsefeats.jl")
include("features/features.jl")

include("radius.jl")
include("transformation.jl")

include("gmms.jl")

include("conformers/bondrotate.jl")
include("conformers/coarsealign.jl")
include("conformers/conformers.jl")

include("draw.jl")

end
