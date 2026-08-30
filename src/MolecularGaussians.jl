"""
MolecularGaussians.jl
===========================

MolecularGaussians.jl is a package used to align molecules and their pharmacophore features by modeling them as Gaussian mixture models.
It makes use of GaussianMixtureAlignment.jl to compute alignments, overlap, and distances between molecules.

REPL help
=========

`?` followed by an algorithm or constructor name prints its help to the terminal, e.g.

- `?PharmacophoreGMM`
- `?gogma_align`
- `?tiv_gogma_align`
- `?rocs_align`
"""
module MolecularGaussians

using StaticArrays: SVector
using CoordinateTransformations: LinearMap, Translation
using Rotations: AngleAxis, RotMatrix
using LinearAlgebra: normalize

using MolecularGraph: MolecularGraph, SDFAtom, SDFMolGraph, SimpleMolGraph, get_prop, is_rotatable, moldisplay, molecular_formula, props, set_prop!, smartstomol, substruct_matches
using Graphs: Graphs, connected_components, edges, induced_subgraph, neighbors, vertices

# MolecularGraph renamed `atomnumber` to `atom_number` in v0.19.0.
@static if pkgversion(MolecularGraph) < v"0.19.0"
    using MolecularGraph: atomnumber as atom_number
else
    using MolecularGraph: atom_number
end

using GaussianMixtureAlignment: IsotropicGaussian, StackedLabeledGaussian, StackedLabeledIsotropicGMM, AbstractGMM, AbstractStackedLabeledIsotropicGMM, ROCSAlignmentResult, centroid, local_align, rocs_align, gogma_align, tiv_gogma_align, overlap, distance, tanimoto
using GaussianMixtureAlignment: flex_gogma_align, aligned, joint_angles
# imported (not just `using`) because the PharmacophoreGMM methods below extend them
import GaussianMixtureAlignment: njoints, joint_axis, joint_origin, joint_features, joint_children, flex

# The alignment functions default to AutoForwardDiff(); loading ForwardDiff
# activates the DifferentiationInterface backend that default requires.
import ForwardDiff

export local_align, gogma_align, tiv_gogma_align, overlap, distance, tanimoto
export PharmacophoreGMM, feature_labels, bondrotate
export flexible_align, flex_gogma_align, aligned, joint_angles
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

Returns the `FigureAxisPlot` produced by the plot, which can be used to further
customize, display, or save the figure.
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
include("conformers/flexiblealign.jl")

end
