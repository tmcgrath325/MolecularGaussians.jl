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

using MolecularGraph: MolecularGraph, SDFAtom, SDFMolGraph, SimpleMolGraph, atom_radius, atomnumber, get_prop, is_rotatable, moldisplay, molecular_formula, props, smartstomol, substruct_matches
using Graphs: Graphs, connected_components, edges, induced_subgraph, neighbors, vertices

using GaussianMixtureAlignment: IsotropicGaussian, IsotropicGMM, AbstractIsotropicMultiGMM, centroid, local_align, rocs_align, gogma_align, tiv_gogma_align, overlap, distance, tanimoto

# The alignment functions default to AutoForwardDiff(); loading ForwardDiff
# activates the DifferentiationInterface backend that default requires.
import ForwardDiff

export local_align, gogma_align, tiv_gogma_align, overlap, distance, tanimoto
export PharmacophoreGMM
export rocs_align

export AtomType, FeatureDef
export parse_feature_definitions, feature_maps

export moldisplay

"""
    pharmacophoredisplay(gmms...; kwargs...)

Plot one or more `PharmacophoreGMM`s and their underlying molecules.

Requires a Makie backend (e.g. GLMakie or CairoMakie) to be loaded; the method
is provided by the MakieCore package extension. Without a backend loaded, this
call raises a `MethodError`.
"""
function pharmacophoredisplay end

# `public` is a soft keyword only on Julia 1.11+; guard so the module still
# loads on the 1.10 LTS, where the declaration is simply absent.
@static if VERSION >= v"1.11"
    eval(Meta.parse("public pharmacophoredisplay"))
end

include("utils.jl")

include("features/featuredef.jl")
include("features/parsefeats.jl")
include("features/features.jl")

include("radius.jl")
include("transformation.jl")

include("gmms.jl")

include("conformers/bondrotate.jl")
include("conformers/conformers.jl")

end
