using MolecularGaussians
using MolecularGraph
using StaticArrays
using CoordinateTransformations
using Rotations
using LinearAlgebra
using Test

using Graphs: induced_subgraph
using GaussianMixtureAlignment: distance
using MolecularGaussians: nodeset

const MG = MolecularGaussians
const FAMILY_DEFS = parse_feature_definitions()


@testset "Gaussian Mixture Distance" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    gonane = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "gonane5α.sdf"))
    remove_hydrogens!(gonane)
    # identical mixture models have no distance
    mol_fmaps = feature_maps(mol, FAMILY_DEFS, [:Volume])
    mol_gmm = PharmacophoreGMM(mol; featuremaps = mol_fmaps)
    gonane_fmaps = feature_maps(gonane, FAMILY_DEFS, [:Volume])
    gonane_gmm = PharmacophoreGMM(gonane; featuremaps = gonane_fmaps)
    @test abs(distance(mol_gmm, mol_gmm)) < 1e-12
    @test abs(distance(gonane_gmm, gonane_gmm)) < 1e-12
    # different mixture models have some distance
    @test distance(mol_gmm, gonane_gmm) > 0.1
    # different transforms of a molecule have some distance
    gmm = PharmacophoreGMM(mol)
    tform = AffineMap(RotationVec(π*rand(3)...), SVector(5*rand(3)...))
    @test distance(tform(gmm), gmm) > 0.1
    # subgraphs of a molecule have some distance
    submol, _ = induced_subgraph(mol, collect(nodeset(mol))[1:Int(floor(end/2))])
    sub_gmm = PharmacophoreGMM(submol)
    @test distance(sub_gmm, gmm) > 0.1
end

@testset "PharmacophoreGMM alignment" begin
    ## PharmacophoreGMM alignment
    mol1 = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    mol2 = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1103_3d.sdf"))
    molgmm1 = PharmacophoreGMM(mol1)
    molgmm2 = PharmacophoreGMM(mol2)
    gonane = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "gonane5α.sdf"))
    # No rotation/translation when aligning a mol to itself
    self_res = rocs_align(molgmm1, molgmm1)
    @test norm(self_res.tform.translation) ≈ 0 atol=1e-6
    @test self_res.tform.linear ≈ I
    # Do you get similar distances when performing the alignment in both directions?
    f_res = rocs_align(molgmm1, molgmm2)
    b_res = rocs_align(molgmm2, molgmm1)
    f_ovrlp, b_ovrlp = f_res.minimum, b_res.minimum
    @test abs(2*(f_ovrlp-b_ovrlp)/(f_ovrlp+b_ovrlp)) < 0.01
end