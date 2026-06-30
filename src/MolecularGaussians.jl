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

using StaticArrays: SVector
using CoordinateTransformations: AffineMap, LinearMap, Translation
using Rotations: AngleAxis, RotMatrix

using MolecularGraph: MolecularGraph, SDFAtom, SDFMolGraph, SimpleMolGraph, atomnumber, colortype, Color, get_prop, is_rotatable, moldisplay, molecular_formula, props, smartstomol, stick!, substruct_matches
using Graphs: Graphs, connected_components, edges, induced_subgraph, neighbors, vertices

using GaussianMixtureAlignment: IsotropicGaussian, IsotropicGMM, AbstractIsotropicMultiGMM, centroid, local_align, rocs_align, gogma_align, tiv_gogma_align, overlap, distance, tanimoto, gmmdisplay!

# The alignment functions default to AutoForwardDiff(); loading ForwardDiff
# activates the DifferentiationInterface backend that default requires.
import ForwardDiff

export local_align, gogma_align, tiv_gogma_align, overlap, distance, tanimoto
export PharmacophoreGMM
export rocs_align

export AtomType, FeatureDef
export parse_feature_definitions, feature_maps

export moldisplay

using Colors: Colors

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
