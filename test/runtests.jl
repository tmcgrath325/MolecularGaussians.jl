using MolecularGaussians
using GaussianMixtureAlignment
using MolecularGraph
using MolecularGraph.Graph
using LinearAlgebra
using StaticArrays
using CoordinateTransformations
using Optim
using Test

using GaussianMixtureAlignment: distance

const MG = MolecularGaussians

@testset "Partial Charge" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "data", "E1050_3d.sdf"))
    # adding partial charges to attributes
    @test_throws KeyError pc = mol.attributes[:partialcharges]
    partialcharges!(mol)
    @test typeof(mol.attributes[:partialcharges]) <: Dict
    # A sulphur atom should have the most positive partial charge
    @test nodeattr(mol,findmax(partialcharges!(mol))[2]).symbol == :S
    # An oxygen atom should have the most negative partial charge
    @test nodeattr(mol,findmin(partialcharges!(mol))[2]).symbol == :O
end

@testset "Van der Waals radii" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "data", "E1050_3d.sdf"))
    # adding radii to attributes
    @test_throws KeyError pc = mol.attributes[:vdwradii]
    vdwradii!(mol)
    @test typeof(mol.attributes[:vdwradii]) <: Dict
    # All atoms should have the appropriate radius
    for idx in nodeset(mol)
        @test vdwradii!(mol)[idx] == MolecularGraph.ATOM_VANDERWAALS_RADII[atomnumber(nodeattr(mol,idx).symbol)]
    end
end

@testset "Gaussian Mixture Distance" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "data", "E1050_3d.sdf"))
    gonane = sdftomol(joinpath(@__DIR__, "..", "data", "gonane5α.sdf"))
    # identical mixture models have no distance
    mol_gmm = MolGMM(mol)
    gonane_gmm = MolGMM(gonane)
    @test abs(distance(mol_gmm, mol_gmm)) < 1e-12
    @test abs(distance(gonane_gmm, gonane_gmm)) < 1e-12
    # different mixture models have some distance
    @test distance(mol_gmm, gonane_gmm) > 0.1
    # different transforms of a molecule have some distance
    gmm= MolGMM(mol)
    tform = AffineMap(10*rand(6)...)
    @test distance(tform(gmm), gmm) > 0.1
    # subgraphs of a molecule have some distance
    submol = nodesubgraph(mol, collect(nodeset(mol))[1:Int(floor(end/2))])
    sub_gmm = MolGMM(submol)
    @test distance(sub_gmm, gmm) > 0.1
end

@testset "MolGMM and PharmacophoreGMM alignment" begin
    ## MolGMM alignment
    mol1= sdftomol(joinpath(@__DIR__, "..", "data", "E1050_3d.sdf"))
    mol2 = sdftomol(joinpath(@__DIR__, "..", "data", "E1103_3d.sdf"))
    molgmm1 = MolGMM(mol1)
    molgmm2 = MolGMM(mol2)
    gonane = sdftomol(joinpath(@__DIR__, "..", "data", "gonane5α.sdf"))
    # No rotation/translation when aligning a mol to itself
    self_res = tiv_gogma_align(molgmm1, molgmm1, 0.5, 0.5; maxstagnant=1000)
    self_tform = AffineMap(self_res.tform_params...)
    @test norm(self_tform.translation) ≈ 0 atol=1e-6
    @test self_tform.linear ≈ I
    # Do you get similar distances when performing the alignment in both directions?
    f_ovrlp = tiv_gogma_align(molgmm1, molgmm2, 0.5, 0.5; maxstagnant=1000).upperbound
    b_ovrlp = tiv_gogma_align(molgmm2, molgmm1, 0.5, 0.5; maxstagnant=1000).upperbound
    @test abs(2*(f_ovrlp-b_ovrlp)/(f_ovrlp+b_ovrlp)) < 0.05     # 5% difference too forgiving?
    # Do you get the inverse transformation when swapping molfixed and molmoving (requires more evals)?
    # This may not be important for generating a distance matrix for many molecules
    # f_res= tiv_gogma_align(PharmacophoreGMM(mol1), PharmacophoreGMM(mol2); maxstagnant=1000)
    # f_tform = AffineMap(f_res.tform_params...)
    # b_res = tiv_gogma_align(PharmacophoreGMM(mol2), PharmacophoreGMM(mol1); maxstagnant=1000)
    # b_tform = AffineMap(b_res.tform_params...)
    # @test (b_tform ∘ f_tform).linear ≈ I atol=0.1      
    # @test (b_tform ∘ f_tform).tranlsation ≈ zeros(3) atol=0.1     
end