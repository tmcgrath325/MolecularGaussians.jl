using MolecularGraph: smilestomol

# Gasteiger charges folded onto heavy atoms, for comparison with RDKit: MolecularGraph makes
# bracket hydrogens ([NH3+]) explicit vertices, while RDKit reports their charge on the heavy
# atom.
function _heavy_group_charges(mol)
    qg = MG.atom_group_charges(mol)
    syms = MG.atom_symbol(mol)
    return [qg[i] for i in eachindex(qg) if syms[i] != :H]
end

@testset "Gasteiger partial charges" begin
    # Reference values from RDKit (ComputeGasteigerCharges, nIter=12, charge + implicit-H
    # charge per heavy atom), the de facto reference implementation.
    rdkit = [
        "CCO" => [0.034281, 0.152360, -0.186641],
        "CN" => [0.097074, -0.097074],
        "CC(=O)O" => [0.138313, 0.299685, -0.252820, -0.185177],
        "CF" => [0.254730, -0.254730],
        "CS(C)(=O)=O" => [0.157362, 0.144104, 0.157362, -0.229414, -0.229414],
        "c1ccccc1Cl" => [0.020312, 0.001534, 0.000106, 0.001534, 0.020312, 0.040548, -0.084344],
        "CC#N" => [0.139944, 0.058715, -0.198658],
        "CS(C)=O" => [0.122677, 0.014804, 0.122677, -0.260158],
        "c1ccsc1" => [0.011626, 0.011626, 0.064602, -0.152454, 0.064602],
        "CP(C)C" => [0.038779, -0.116337, 0.038779, 0.038779],
        "CC(=O)[O-]" => [0.062675, 0.038279, -0.550477, -0.550477],   # split across the carboxylate oxygens
        "C[NH3+]" => [0.328678, 0.671322],                            # localized ammonium charge
        "C[N+](=O)[O-]" => [0.513391, 0.016283, -0.264837, -0.264837],
    ]
    for (smi, ref) in rdkit
        q = _heavy_group_charges(smilestomol(smi))
        @test q ≈ ref atol = 1e-5
        @test sum(q) ≈ round(sum(ref)) atol = 1e-10   # charges sum to the total formal charge
    end

    # Sulfonate: the seed charge splits equally across the three equivalent oxygens (RDKit
    # leaves it on the charged one), so the two sp2 oxygens agree and the total is exact.
    qs = _heavy_group_charges(smilestomol("CS(=O)(=O)[O-]"))
    @test qs[3] ≈ qs[4]
    @test sum(qs) ≈ -1 atol = 1e-10

    # explicit hydrogens (SDF input) are charged like any atom and folded by group charges
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    q = partial_charges(mol)
    syms = MG.atom_symbol(mol)
    @test length(q) == length(vertices(mol))
    @test sum(q) ≈ 0 atol = 1e-10
    @test all(q[i] > 0 for i in eachindex(q) if syms[i] == :H)   # H bound to C/O is positive here
    qg = MG.atom_group_charges(mol)
    @test all(iszero(qg[i]) for i in eachindex(qg) if syms[i] == :H)
    @test sum(qg) ≈ 0 atol = 1e-10
    # sulfonate ester sulfurs carry the most positive group charge in this molecule
    @test all(qg[i] > 0.3 for i in eachindex(qg) if syms[i] == :S)

    # elements outside the parameterization fail fast
    @test_throws "no Gasteiger parameters" partial_charges(smilestomol("C[SeH]"))
end

@testset "physics fields" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    heavy = MG.heavy_atom_idxs(mol)
    qg = MG.atom_group_charges(mol)

    p = physics_gmm(mol)
    @test p isa PharmacophoreGMM{3,Float64,3,Symbol}
    @test length(p) == length(heavy)
    @test sort(feature_labels(p)) == sort([:Charge, :Steric, :VdW])
    @test GMA.njoints(p) > 0                       # rotatable bonds recorded as usual
    @test GMA.njoints(physics_gmm(mol; rigid=true)) == 0

    # slot values follow the per-family conventions
    for g in p.gaussians
        i = heavy[findfirst(h -> props(mol, h).coords == g.μ, heavy)]
        atom = props(mol, i)
        slot = Dict(l => k for (k, l) in enumerate(g.labels))
        @test g.σ[slot[:VdW]] ≈ MG.vdw_volume_sigma(atom, MG.rocs_amplitude) && g.ϕ[slot[:VdW]] ≈ MG.rocs_amplitude
        @test g.σ[slot[:Steric]] ≈ 0.5 * MG.vdw_radius(atom) && g.ϕ[slot[:Steric]] ≈ MG.rocs_amplitude
        @test g.σ[slot[:Charge]] ≈ 1.0 && g.ϕ[slot[:Charge]] ≈ qg[i]
    end
    # parameters flow through
    p2 = physics_gmm(mol; steric_scale=0.3, charge_sigma=2.0, amplitude=1.0)
    g2 = p2.gaussians[1]
    slot2 = Dict(l => k for (k, l) in enumerate(g2.labels))
    @test g2.ϕ[slot2[:VdW]] == 1.0 && g2.σ[slot2[:Charge]] == 2.0

    # interaction coefficients: attraction for :VdW and :Charge, repulsion for :Steric
    inter = physics_interactions(; steric=0.5)
    @test inter[(:VdW, :VdW)] == 1.0 && inter[(:Steric, :Steric)] == -0.5 && inter[(:Charge, :Charge)] == 1.0

    # the fields compose: overlap and rigid local alignment run with the physics coefficients
    q = physics_gmm(sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1103_3d.sdf")))
    o = overlap(p, q; interactions=inter)
    @test isfinite(o)
    # matching-sign charge overlap beats the sign-flipped one on identical copies
    pσ, pϕ = GMA.pairwise_consts(p, p, physics_interactions())
    pσn, pϕn = GMA.pairwise_consts(p, p, physics_interactions(; charge=-1.0))
    @test overlap(p, p, pσ, pϕ) > overlap(p, p, pσn, pϕn)

    # a families subset skips the charge computation and the absent labels
    pv = physics_gmm(mol; families=(:VdW, :Steric))
    @test sort(feature_labels(pv)) == [:Steric, :VdW]
end
