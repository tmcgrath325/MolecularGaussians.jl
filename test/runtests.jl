using MolecularGaussians
using MolecularGraph
using MolecularGraph.Graph
using LinearAlgebra
using StaticArrays
using CoordinateTransformations
using Optim
using Test

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
    @test abs(gmm_distance(mol_gmm, mol_gmm)) < 1e-12
    @test abs(gmm_distance(gonane_gmm, gonane_gmm)) < 1e-12
    # different mixture models have some distance
    @test gmm_distance(mol_gmm, gonane_gmm) > 0.1
    # different transforms of a molecule have some distance
    gmmfixed = MolGMM(mol)
    gmmmoving = MolGMM(mol)
    tform = AffineMap( [0. 1. 0.;
                        1. 0. 0.;
                        0. 0. 1.],
                       [1., 0., 0.]
    )
    @test gmm_distance(gmmfixed, gmmmoving; tform2=tform) > 0.1
    # subgraphs of a molecule have some distance
    molscs = nodesubgraph(mol, collect(nodeset(mol))[1:Int(floor(end/2))])
    molscs_gmm = MolGMM(molscs)
    @test gmm_distance(mol_gmm, molscs_gmm) > 0.1
end

@testset "GMM alignment" begin
    ## MolGMM alignment
    mol1= sdftomol(joinpath(@__DIR__, "..", "data", "E1050_3d.sdf"))
    mol2 = sdftomol(joinpath(@__DIR__, "..", "data", "E1103_3d.sdf"))
    molgmm1 = MolGMM(mol1)
    molgmm2 = MolGMM(mol2)
    gonane = sdftomol(joinpath(@__DIR__, "..", "data", "gonane5α.sdf"))
    # No rotation/translation when aligning a mol to itself
    selftform = tiv_gogma_align_gmms(molgmm1, molgmm1, 2, 2; maxstagnant=1000)[2]
    @test norm(selftform.translation) ≈ 0 atol=1e-6
    @test selftform.linear ≈ I
    # Do you get similar distances when performing the alignment in both directions?
    fovrlp, frwdtform = tiv_gogma_align_gmms(molgmm1, molgmm2, 2, 2; maxstagnant=1000)[1:2]
    bovrlp, bkwdtform = tiv_gogma_align_gmms(molgmm2, molgmm1, 2, 2; maxstagnant=1000)[1:2]
    @test abs(2*(fovrlp-bovrlp)/(fovrlp+bovrlp)) < 0.05     # 5% difference too forgiving?
    # Do you get the inverse transformation when swapping molfixed and molmoving (requires more evals)?
    # This may not be important for generating a distance matrix for many molecules
    fovrlp, frwdtform = tiv_gogma_align_gmms(FeatureMolGMM(mol1), FeatureMolGMM(mol2))[1:2]
    bovrlp, bkwdtform = tiv_gogma_align_gmms(FeatureMolGMM(mol2), FeatureMolGMM(mol1))[1:2]
    @test frwdtform.linear * bkwdtform.linear ≈ I atol=0.1      
end