using MolecularGaussians
using MolecularGraph
using StaticArrays
using CoordinateTransformations
using Rotations
using Test

using Graphs: induced_subgraph
using GaussianMixtureAlignment: distance
using MolecularGaussians: nodeset

const MG = MolecularGaussians

@testset "Gaussian Mixture Distance" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    gonane = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "gonane5α.sdf"))
    # identical mixture models have no distance
    mol_gmm = MolGMM(mol)
    gonane_gmm = MolGMM(gonane)
    @test abs(distance(mol_gmm, mol_gmm)) < 1e-12
    @test abs(distance(gonane_gmm, gonane_gmm)) < 1e-12
    # different mixture models have some distance
    @test distance(mol_gmm, gonane_gmm) > 0.1
    # different transforms of a molecule have some distance
    gmm= MolGMM(mol)
    tform = AffineMap(RotationVec(π*rand(3)...), SVector(5*rand(3)...))
    @test distance(tform(gmm), gmm) > 0.1
    # subgraphs of a molecule have some distance
    submol, _ = induced_subgraph(mol, collect(nodeset(mol))[1:Int(floor(end/2))])
    sub_gmm = MolGMM(submol)
    @test distance(sub_gmm, gmm) > 0.1
end

@testset "MolGMM and PharmacophoreGMM alignment" begin
    ## MolGMM alignment
    mol1 = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    mol2 = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1103_3d.sdf"))
    molgmm1 = MolGMM(mol1)
    molgmm2 = MolGMM(mol2)
    gonane = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "gonane5α.sdf"))
    # No rotation/translation when aligning a mol to itself
    self_res = rocs_align(molgmm1, molgmm1)
    @test norm(self_res.tform.translation) ≈ 0 atol=1e-6
    @test self_res.tform.linear ≈ I
    # Do you get similar distances when performing the alignment in both directions?
    f_res = rocs_align(molgmm1, molgmm2)
    b_res = rocs_align(molgmm2, molgmm1)
    f_ovrlp, b_ovrlp = f_res.upperbound, b_res.upperbound
    @test abs(2*(f_ovrlp-b_ovrlp)/(f_ovrlp+b_ovrlp)) < 0.01
end