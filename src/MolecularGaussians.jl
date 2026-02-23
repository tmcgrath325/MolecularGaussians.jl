"""
MolecularGaussians.jl
===========================

MolecularGaussians.jl is a package used to align molecules and their pharmacophore features by modeling them as Gaussian mixture models.
It makes use of GaussianMixtureAlignment.jl to compute alignments, overlap, and distances between molecules.

REPL help
=========

? followed by an algorith or constructor name will print help to the terminal. See: \n
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
using MolecularGraph: SimpleMolGraph
using Graphs

using Optim

using GaussianMixtureAlignment
using GaussianMixtureAlignment: AbstractGaussian, AbstractSingleGMM, AbstractMultiGMM, AbstractGMM
using GaussianMixtureAlignment: AbstractIsotropicGaussian, AbstractIsotropicGMM, AbstractIsotropicMultiGMM
using GaussianMixtureAlignment: IsotropicGaussian, IsotropicGMM, IsotropicMultiGMM
using GaussianMixtureAlignment: centroid
using GaussianMixtureAlignment: local_align, rocs_align, gogma_align, rot_gogma_align, tiv_gogma_align, overlap, distance, tanimoto

export PharmacophoreGMM

export AtomType, FeatureDef
export parse_feature_definitions, feature_maps

export rotatablesubgraphs, rotatablebonds, conformers
export rotate_edge, rotate_edge!, rotate_edges, rotate_edges!
export nodeset

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
include("conformers/conformers.jl")

include("draw.jl")

end
