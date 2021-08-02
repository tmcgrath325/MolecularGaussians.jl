module MolecularGaussians

using MolecularGraph
using MolecularGraph.Graph
using MolecularGraph: atommatch, bondmatch, emaptonmap
using CoordinateTransformations
using StaticArrays

using LinearAlgebra
using Optim
using GOGMA

using PlotlyJS

export MolGMM, FeatureMolGMM, gmm_overlap, gmm_distance
export vdwradii, vdwradii!
export partialcharges, partialcharges!
export pharmfeatures, pharmfeatures!
export local_align_gmms
export inertial_transforms, rocs_align_gmms
export gogmaGMM, gogma_align_gmms, tiv_gogma_align_gmms
export drawMolGMM, drawMolGMMs, drawmol, drawFeatureMolGMMs, plotdrawing

include("utils.jl")
include("gmms/gmms.jl")
include("radius.jl")
include("partialcharge.jl")
include("pharmfeatures.jl")
include("gmms/pharmacophores.jl")
include("gmms/distance.jl")
include("gmms/align.jl")
include("gmms/rocsalign.jl")
include("gmms/gogma.jl")
include("draw.jl")
end
