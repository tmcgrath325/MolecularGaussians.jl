import GaussianMixtureAlignment as GMA
import ForwardDiff

@testset "flexible alignment" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1103_3d.sdf"))
    p = PharmacophoreGMM(mol)
    K = GMA.njoints(p)
    @test K == length(p.axes) == 2

    # the articulated-model interface reflects the stored rotatable-bond geometry
    for b in 1:K
        @test GMA.joint_axis(p, b) ≈ p.axes[b] / norm(p.axes[b])   # unit direction
        @test norm(GMA.joint_axis(p, b)) ≈ 1
        @test GMA.joint_origin(p, b) == p.origins[b]
        @test GMA.joint_features(p, b) == p.bondtogaussians[b]
        @test GMA.joint_children(p, b) == p.bondtobonds[b]
    end
    # each rotatable bond precedes the bonds downstream of it — the ordering the search assumes
    @test all(all(c > b for c in GMA.joint_children(p, b)) for b in 1:K)

    # flex is the forward kinematics `bondrotate` already provides, and is the identity at zero
    φ = [0.3, -0.5]
    @test all(GMA.flex(p, φ).gaussians[i].μ ≈ bondrotate(p, φ, 1:K).gaussians[i].μ for i in eachindex(p.gaussians))
    @test all(GMA.flex(p, zeros(K)).gaussians[i].μ ≈ p.gaussians[i].μ for i in eachindex(p.gaussians))

    # bond rotation carries ForwardDiff duals, so the alignment objective is differentiable
    obj(φ) = sum(abs2, GMA.flex(p, φ).gaussians[3].μ)
    @test all(isfinite, ForwardDiff.gradient(obj, [0.2, -0.3]))

    # a continuous flexible search does at least as well as brute-force enumeration of a coarse
    # conformer grid aligned rigidly — the payoff over the discrete `conformers` workflow
    target = GMA.flex_pose((0.3, -0.2, 0.4, 0.8, -1.0, 0.5, 0.9, -0.7), p)
    confs = PharmacophoreGMM.(MG.conformers(mol; step = Float64(π), maxbonds = 2))
    _, _, _, brute_obj = MG.align_conformers(confs, target; alignfun = local_align)
    res = flexible_align(p, target; maxsplits = 40)
    @test -res.upperbound >= -brute_obj - 1.0e-6
    @test res.lowerbound <= res.upperbound

    # the result carries the posed, flexed molecule and one angle per rotatable bond
    @test aligned(res) isa PharmacophoreGMM
    @test length(joint_angles(res)) == K
end

@testset "flexible alignment: articulated target" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1103_3d.sdf"))
    p = PharmacophoreGMM(mol)
    K = GMA.njoints(p)
    target = GMA.flex_pose((0.3, -0.2, 0.4, 0.8, -1.0, 0.5, 0.9, -0.7), p)

    # a molecule used as the target keeps its input conformation unless asked otherwise
    fixed = flexible_align(p, target; maxsplits = 10)
    @test GMA.target_joint_angles(fixed) == ()
    @test all(GMA.aligned_target(fixed).gaussians[i].μ == target.gaussians[i].μ for i in eachindex(target.gaussians))

    # with `flextarget`, the target's rotatable bonds are searched too
    both = flexible_align(p, target; flextarget = true, maxsplits = 10)
    @test length(GMA.target_joint_angles(both)) == K
    @test GMA.aligned_target(both) isa PharmacophoreGMM
    @test both.lowerbound <= both.upperbound

    # the self-overlap penalty leaves the reported objective consistent with the posed models
    pen = flexible_align(p, target; selfoverlap = 1.0, maxsplits = 10)
    xt = aligned(pen)
    @test pen.upperbound ≈ -overlap(xt, target) + GMA.penalty(GMA.SelfOverlap(p), xt) atol = 1.0e-8
end
