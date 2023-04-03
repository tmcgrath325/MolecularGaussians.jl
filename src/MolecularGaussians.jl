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
using MolecularGraph.Graph
using MolecularGraph: atommatch, bondmatch, emaptonmap

using GaussianMixtureAlignment
using GaussianMixtureAlignment: AbstractGaussian, AbstractSingleGMM, AbstractMultiGMM, AbstractGMM
using GaussianMixtureAlignment: AbstractIsotropicGaussian, AbstractIsotropicGMM, AbstractIsotropicMultiGMM
using GaussianMixtureAlignment: IsotropicGaussian, IsotropicGMM, IsotropicMultiGMM
using GaussianMixtureAlignment: centroid
using GaussianMixtureAlignment: local_align, rocs_align, gogma_align, rot_gogma_align, tiv_gogma_align, overlap, distance, tanimoto

export local_align, gogma_align, tiv_gogma_align, overlap, distance, tanimoto
export add_attributes!, attributes
export AtomGaussian, MolGMM, PharmacophoreGMM
export vdwradii, vdwradii!
export partialcharges, partialcharges!
export pharmfeatures, pharmfeatures!
export inertial_transforms, rocs_align
export affinetransform

using MakieCore
using GeometryBasics
using Colors

include("utils.jl")
include("gmms/gmms.jl")
include("transformation.jl")
include("radius.jl")
include("partialcharge.jl")
include("pharmfeatures.jl")
include("gmms/pharmacophores.jl")

include("draw/drawmol.jl")

end
