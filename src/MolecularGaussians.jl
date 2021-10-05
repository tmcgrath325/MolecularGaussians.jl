module MolecularGaussians

using LinearAlgebra

using StaticArrays
using CoordinateTransformations

using MolecularGraph
using MolecularGraph.Graph
using MolecularGraph: atommatch, bondmatch, emaptonmap

using GaussianMixtureAlignment
using GaussianMixtureAlignment: AbstractGaussian, AbstractSingleGMM, AbstractMultiGMM, AbstractGMM
using GaussianMixtureAlignment: AbstractIsotropicGaussian, AbstractIsotropicGMM, AbstractIsotropicMultiGMM
using GaussianMixtureAlignment: IsotropicGaussian, IsotropicGMM, IsotropicMultiGMM
using GaussianMixtureAlignment: center_of_mass
using GaussianMixtureAlignment: local_align, gogma_align, rot_gogma_align, tiv_gogma_align, overlap, distance, tanimoto

using PlotlyJS

export local_align, gogma_align, tiv_gogma_align, overlap, distance, tanimoto
export add_attributes!, attributes
export AtomGaussian, MolGMM, PharmacophoreGMM
export vdwradii, vdwradii!
export partialcharges, partialcharges!
export pharmfeatures, pharmfeatures!
export inertial_transforms, rocs_align
export drawMolGMM, drawMolGMMs, drawmol, drawPharmacophoreGMMs, plotdrawing

include("utils.jl")
include("gmms/gmms.jl")
include("transformation.jl")
include("radius.jl")
include("partialcharge.jl")
include("pharmfeatures.jl")
include("gmms/pharmacophores.jl")
include("gmms/rocsalign.jl")
include("draw.jl")
end
