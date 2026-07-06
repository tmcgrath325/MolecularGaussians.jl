using MolecularGaussians
using MolecularGraph
using StaticArrays
using CoordinateTransformations
using Rotations
using LinearAlgebra
using Test
using Aqua
using ExplicitImports
using MakieCore   # triggers the MolecularGaussiansMakieCoreExt package extension
using OffsetArrays: OffsetArray
using Documenter

using Graphs: induced_subgraph, edges, vertices, connected_components, rem_edge!, neighbors
using GaussianMixtureAlignment: distance, IsotropicGaussian, LabeledIsotropicGMM
using MolecularGaussians: nodeset

const MG = MolecularGaussians
const FAMILY_DEFS = parse_feature_definitions()

DocMeta.setdocmeta!(MolecularGaussians, :DocTestSetup,
                    :(using MolecularGaussians, MolecularGraph); recursive = true)

@testset "Doctests" begin
    doctest(MolecularGaussians)
end

include("flexible.jl")


@testset "Gaussian Mixture Distance" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    gonane = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "gonane5α.sdf"))
    remove_hydrogens!(gonane)
    # identical mixture models have no distance
    mol_fmaps = feature_maps(mol, FAMILY_DEFS, [:Volume])
    mol_gmm = PharmacophoreGMM(mol; feature_maps = mol_fmaps)
    gonane_fmaps = feature_maps(gonane, FAMILY_DEFS, [:Volume])
    gonane_gmm = PharmacophoreGMM(gonane; feature_maps = gonane_fmaps)
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

@testset "PharmacophoreGMM keyword constructor" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    # σfun and ϕfun are keyword arguments; the model stores what was passed
    pgmm = PharmacophoreGMM(mol; ϕfun = a -> 2.0)
    @test pgmm isa PharmacophoreGMM
    @test pgmm.ϕfun(props(mol, 1)) == 2.0
    # the former positional σfun/ϕfun form is a clean break — no shim
    @test_throws MethodError PharmacophoreGMM(mol, a -> 2.0)
    # `feature_maps` is both the keyword and the field name (matching the function)
    fm = feature_maps(mol, FAMILY_DEFS, [:Volume])
    pgmm_fm = PharmacophoreGMM(mol; feature_maps = fm)
    @test pgmm_fm.feature_maps == fm
    # the type-based eltype agrees with the instance-based one: a PharmacophoreGMM
    # is a flat stacked GMM, so its element type is the Gaussian, not a keyed pair
    @test eltype(typeof(pgmm)) == eltype(pgmm) == MG.StackedLabeledGaussian{3,Float64,1,Symbol}
end

@testset "labeled PharmacophoreGMM" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    pgmm = PharmacophoreGMM(mol)
    # one component per feature when no two features share a set of atoms
    @test length(pgmm.gaussians) == length(pgmm) == 28
    @test feature_labels(pgmm) == [:Volume]
    # self-overlap is a regression pin: identity-label interactions reproduce the
    # per-family "only same-key GMMs score" overlap.
    @test MG.overlap(pgmm, pgmm) ≈ 171.00457064028473 rtol=1e-10

    # cross-label interactions are a new capability the per-family Dict design
    # could not express: distinct families only overlap when paired explicitly.
    fm = feature_maps(mol, FAMILY_DEFS, [:Donor, :Acceptor])
    pg = PharmacophoreGMM(mol; feature_maps = fm)
    same = MG.overlap(pg, pg)                                        # equal labels only
    crossed = MG.overlap(pg, pg; interactions = Dict((:Donor, :Acceptor) => 1.0))
    @test same != crossed
end

@testset "stacked PharmacophoreGMM" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    families = [:Donor, :Acceptor, :Aromatic, :NegIonizable]
    fm = feature_maps(mol, FAMILY_DEFS, families)
    pgmm = PharmacophoreGMM(mol; feature_maps = fm)
    nfeatures = sum(length, values(fm))

    # Features built from the same atoms — here hydroxyls, which are both donors and
    # acceptors — share one component, so there are fewer components than features and
    # the stacking degree is the largest number of features on any one atom set.
    @test length(pgmm.gaussians) < nfeatures
    @test pgmm isa PharmacophoreGMM{3,Float64,2,Symbol}
    @test sum(g -> count(!iszero, g.ϕ), pgmm.gaussians) == nfeatures
    @test sort(feature_labels(pgmm)) == sort(families)
    stacked = filter(g -> count(!iszero, g.ϕ) > 1, pgmm.gaussians)
    @test !isempty(stacked)
    @test all(sort(collect(g.labels)) == [:Acceptor, :Donor] for g in stacked)

    # Stacking is a representation change, not a modeling change: overlaps match those of
    # the mean-duplicated labeled model, with and without interaction weights. The two
    # representations sum the same terms in a different order, so they agree to rounding.
    unstacked = LabeledIsotropicGMM(
        [IsotropicGaussian(g.μ, g.σ[j], g.ϕ[j]) for g in pgmm.gaussians
             for j in eachindex(g.labels) if !iszero(g.ϕ[j])],
        [l for g in pgmm.gaussians for (j, l) in pairs(g.labels) if !iszero(g.ϕ[j])])
    @test length(unstacked.gaussians) == nfeatures
    @test MG.overlap(pgmm, pgmm) ≈ MG.overlap(unstacked, unstacked) rtol=1e-14
    interactions = Dict((:Donor, :Acceptor) => 0.5, (:Donor, :Donor) => 1.0)
    @test MG.overlap(pgmm, pgmm; interactions) ≈ MG.overlap(unstacked, unstacked; interactions) rtol=1e-14

    # Padded slots carry zero amplitude and a positive width, and repeat a label the
    # component already has, so they add no overlap and no new feature label.
    padded = filter(g -> count(!iszero, g.ϕ) < 2, pgmm.gaussians)
    @test !isempty(padded)
    @test all(all(>(0), g.σ) && g.labels[1] == g.labels[2] for g in padded)

    # A Dict of per-family functions sizes and weights each family separately, which is
    # what distinguishes the stacked slots of one component from each other.
    σfun = Dict(f => ((a, ϕ) -> 1.0 + i) for (i, f) in enumerate(families))
    ϕfun = Dict(f => (a -> 1.0 / i) for (i, f) in enumerate(families))
    perfamily = PharmacophoreGMM(mol; feature_maps = fm, σfun, ϕfun)
    g = first(x for x in perfamily.gaussians if count(!iszero, x.ϕ) > 1)
    @test allunique(g.σ) && allunique(g.ϕ)
    # a family the Dict does not cover is an error, not a silent default
    @test_throws KeyError PharmacophoreGMM(mol; feature_maps = fm, ϕfun = Dict(:Donor => a -> 1.0))
end

@testset "PharmacophoreGMM bond rotation" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    pgmm = PharmacophoreGMM(mol)
    sgs = MG.rotablesubgraphs(mol)
    # one bond frame per rotatable bond, and the four bond-geometry vectors agree
    @test length(pgmm.axes) == length(sgs)
    @test length(pgmm.axes) == length(pgmm.origins) == length(pgmm.bondtogaussians) == length(pgmm.bondtobonds)
    # the bond-geometry length invariant is enforced at construction
    @test_throws DimensionMismatch MG.PharmacophoreGMM(pgmm.gaussians,
        pgmm.graph, pgmm.σfun, pgmm.ϕfun, pgmm.feature_maps,
        pgmm.axes, pgmm.origins[1:end-1], pgmm.bondtogaussians, pgmm.bondtobonds)

    # rigid=true records no rotatable bonds
    rig = PharmacophoreGMM(mol; rigid = true)
    @test isempty(rig.axes) && isempty(rig.bondtogaussians)

    # a disconnected subgraph still builds: a rotatable bond splits only its own
    # component (this molecule's first-half induced subgraph is disconnected)
    submol, _ = induced_subgraph(mol, collect(nodeset(mol))[1:Int(floor(end/2))])
    @test length(connected_components(submol)) > 1
    @test PharmacophoreGMM(submol) isa PharmacophoreGMM

    b = 1
    rotated = bondrotate(pgmm, 0.7, b)
    # exactly the Gaussians the bond maps to are moved; the rest are untouched
    moved = [i for i in eachindex(pgmm.gaussians) if !(rotated.gaussians[i].μ ≈ pgmm.gaussians[i].μ)]
    @test sort(moved) ⊆ sort(pgmm.bondtogaussians[b])
    @test all(rotated.gaussians[i].μ == pgmm.gaussians[i].μ for i in eachindex(pgmm.gaussians) if i ∉ pgmm.bondtogaussians[b])
    # self-overlap is maximal, so a bond rotation moves the model away from itself
    @test MG.distance(rotated, pgmm) > 1e-6
    # +θ then -θ recovers the original; a full turn is the identity
    @test MG.distance(bondrotate(rotated, -0.7, b), pgmm) < 1e-10
    @test MG.distance(bondrotate(pgmm, 2π, b), pgmm) < 1e-8

    # rotating the GMM about a bond matches rebuilding it from the graph rotated
    # about that bond's edge
    θ = 0.9
    rebuilt = PharmacophoreGMM(MG.rotate_edge(mol, sgs[b].edge, θ))
    gmmrot = bondrotate(pgmm, θ, b)
    @test all(gmmrot.gaussians[k].μ ≈ rebuilt.gaussians[k].μ for k in eachindex(pgmm.gaussians))

    # the sequence form applies rotations in order; mismatched lengths are rejected
    @test bondrotate(pgmm, [0.3, 0.4], [1, 2]) isa PharmacophoreGMM
    @test_throws DimensionMismatch bondrotate(pgmm, [0.3], [1, 2])

    # combineatoms=false gives one Gaussian per atom, so a multi-atom feature set
    # yields more Gaussians than the combined (default) form
    fm = feature_maps(mol, FAMILY_DEFS, [:Aromatic])
    natoms = sum(sum(length, sets; init = 0) for sets in values(fm); init = 0)
    @test any(length(set) > 1 for set in first(values(fm)))       # there is a multi-atom set
    @test length(PharmacophoreGMM(mol; combineatoms = false, feature_maps = fm).gaussians) == natoms
    @test length(PharmacophoreGMM(mol; feature_maps = fm).gaussians) < natoms
end

@testset "atoms_to_feature defaults" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    nodes = collect(nodeset(mol))
    # The documented defaults (ϕfun=rocs_volume_amplitude, σfun=vdw_volume_sigma)
    # build an IsotropicGaussian for both a single-atom and a multi-atom nodeset.
    @test MG.atoms_to_feature(mol, nodes[1:1]) isa MG.IsotropicGaussian
    @test MG.atoms_to_feature(mol, nodes[1:3]) isa MG.IsotropicGaussian
    # A combined feature is a property of its set of atoms, so listing them in another
    # order must give the identical width — the van der Waals radii come from a Float32
    # table, where accumulating the combined volume would otherwise be order-dependent.
    fwd = MG.atoms_to_feature(mol, nodes[1:5])
    rev = MG.atoms_to_feature(mol, reverse(nodes[1:5]))
    @test fwd.σ == rev.σ && fwd.ϕ == rev.ϕ
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

@testset "align_conformers result dispatch" begin
    mol1 = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf")); remove_hydrogens!(mol1)
    mol2 = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1103_3d.sdf")); remove_hydrogens!(mol2)
    x = PharmacophoreGMM(mol1); y = PharmacophoreGMM(mol2)
    confs = [x]
    # align_conformers normalizes each backend's result by dispatching on the
    # returned value, so the named backends and a wrapper returning one of their
    # results all work. (gogma_align is omitted: the global optimizer is too slow.)
    for alignfun in (local_align, rocs_align, (a, b) -> rocs_align(a, b))
        conf, tform, idx, olap = MG.align_conformers(confs, y; alignfun)
        @test conf === x
        @test idx == 1
        @test olap isa Real && isfinite(olap)
        @test tform !== identity
    end
    # extra keyword arguments are forwarded to alignfun, not dropped — through
    # both the (confs, template) and the (xconfs, yconfs) methods.
    received = Ref{Any}(:unset)
    stub(a, b; myopt = nothing) = (received[] = myopt; (0.0, ntuple(_ -> 0.0, 6)))
    MG.align_conformers(confs, y; alignfun = stub, myopt = 42)
    @test received[] == 42
    received[] = :unset
    MG.align_conformers(confs, [y]; alignfun = stub, myopt = 7)
    @test received[] == 7
    # gogma_align / tiv_gogma_align report their overlap and transform in
    # `upperbound` / `tform_params`; the fallback _overlap_tform normalizes that
    # shape (checked directly to avoid running the slow global optimizer).
    @test MG._overlap_tform((upperbound = 3.5, tform_params = ntuple(zero, 6))) == (3.5, ntuple(zero, 6))
end

@testset "coordinate transforms" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    orig1 = copy(props(mol, 1).coords)
    tform = AffineMap(RotationVec(0.3, -0.2, 0.5), SVector(1.0, 2.0, 3.0))
    # `transformed` maps every atom's coordinates through the tform, preserves the
    # coordinate container type, and leaves the source molecule unmutated.
    mol2 = MG.transformed(tform, mol)
    @test all(props(mol2, i).coords ≈ tform(props(mol, i).coords) for i in vertices(mol))
    @test typeof(props(mol2, 1).coords) == typeof(props(mol, 1).coords)
    @test props(mol, 1).coords == orig1
    # Transforming a PharmacophoreGMM carries the transform into its molecular graph.
    pgmm = PharmacophoreGMM(mol)
    pgmm2 = MG.transform(pgmm, tform)
    @test all(props(pgmm2.graph, i).coords ≈ tform(props(mol, i).coords) for i in vertices(mol))
    # `transform` accepts non-affine transforms too (rotation-only, translation-only).
    lin = LinearMap(RotationVec(0.3, -0.2, 0.5))
    tr = Translation(SVector(1.0, 2.0, 3.0))
    pgmm_lin = MG.transform(pgmm, lin)
    pgmm_tr = MG.transform(pgmm, tr)
    @test all(props(pgmm_lin.graph, i).coords ≈ lin(props(mol, i).coords) for i in vertices(mol))
    @test all(props(pgmm_tr.graph, i).coords ≈ tr(props(mol, i).coords) for i in vertices(mol))
end

@testset "PharmacophoreGMM arithmetic" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    pgmm = PharmacophoreGMM(mol)
    t = SVector(1.0, 2.0, 3.0)
    # Subtracting a translation is adding its negation.
    @test distance(pgmm - t, pgmm + (-t)) < 1e-12
    # Translating away from and back to the origin recovers the original model.
    @test distance((pgmm + t) - t, pgmm) < 1e-12
end

@testset "rotate_edge" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    orig = Dict(i => copy(props(mol, i).coords) for i in vertices(mol))
    # A rotatable bond: removing it splits the molecule into two components. Ring
    # bonds do not, and rotate_edge throws an ArgumentError for them.
    splits_in_two(e) = (g = deepcopy(mol); rem_edge!(g, e); length(connected_components(g)) == 2)
    edge = first(Iterators.filter(splits_in_two, edges(mol)))
    # rotate_edge (non-mutating) leaves the source graph untouched
    rotated = MG.rotate_edge(mol, edge, 0.7)
    @test all(props(mol, i).coords == orig[i] for i in vertices(mol))
    # the smaller side of the bond moves; a full turn returns to the start
    @test any(!(props(rotated, i).coords ≈ orig[i]) for i in vertices(mol))
    @test all(props(MG.rotate_edge(mol, edge, 2π), i).coords ≈ orig[i] for i in vertices(mol))
    # rotate_edge! mutates the graph in place through the public set_prop! path
    g = deepcopy(mol)
    MG.rotate_edge!(g, edge, 0.7)
    @test any(!(props(g, i).coords ≈ orig[i]) for i in vertices(mol))
    # A ring bond does not split the graph in two; rotating about it is undefined.
    ringbond = first(Iterators.filter(!splits_in_two, edges(mol)))
    @test_throws "does not generate two disjoint subgraphs" MG.rotate_edge(mol, ringbond, 0.7)
end

@testset "conformer generation" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    orig = Dict(i => copy(props(mol, i).coords) for i in vertices(mol))
    # rotatablebonds finds the rotatable, non-terminal bonds; rotablesubgraphs
    # wraps each as the RotatableSubgraph that would move if it were rotated.
    rbonds = MG.rotatablebonds(mol)
    @test !isempty(rbonds)
    @test length(MG.rotablesubgraphs(mol)) == length(rbonds)

    # With ignoreH=false, a node is terminal iff it has exactly one neighbor.
    term = first(i for i in vertices(mol) if length(neighbors(mol, i)) == 1)
    inner = first(i for i in vertices(mol) if length(neighbors(mol, i)) > 1)
    @test MG.isterminalnode(mol, term, false)
    @test !MG.isterminalnode(mol, inner, false)

    # conformers enumerates the angle grid over (at most maxbonds) rotatable bonds:
    # `length(lower:step:upper) ^ nbonds` structures. With step=π over [-π, π) that
    # is 2 angles per bond, so a single bond yields 2 conformers.
    confs = MG.conformers(mol; step = π, maxbonds = 1)
    @test length(confs) == 2
    @test all(c isa MolecularGraph.SDFMolGraph for c in confs)
    @test any(any(!(props(c, i).coords ≈ orig[i]) for i in vertices(mol)) for c in confs)

    # rotate_edges applies a sequence of per-bond rotations; the bang form mutates,
    # the plain form leaves the source untouched.
    twobonds = rbonds[1:2]
    rotated = MG.rotate_edges(mol, twobonds, [0.5, 0.7])
    @test all(props(mol, i).coords == orig[i] for i in vertices(mol))
    @test any(!(props(rotated, i).coords ≈ orig[i]) for i in vertices(mol))
    g = deepcopy(mol)
    MG.rotate_edges!(g, twobonds, [0.5, 0.7])
    @test all(props(g, i).coords ≈ props(rotated, i).coords for i in vertices(mol))

    # angleaxis_rotate_coords/_graph rotate about an explicit axis through an origin;
    # a full turn about the z-axis returns every atom to its start.
    zaxis = SVector(0.0, 0.0, 1.0)
    turned = MG.angleaxis_rotate_graph(mol, 2π, zaxis, SVector(0.0, 0.0, 0.0))
    @test turned isa MolecularGraph.SDFMolGraph
    @test all(props(turned, i).coords ≈ orig[i] for i in vertices(mol))
    atom = MG.angleaxis_rotate_coords(2π, zaxis, SVector(0.0, 0.0, 0.0), props(mol, 1))
    @test atom.coords ≈ props(mol, 1).coords
end

@testset "SMARTS atom-type combinators" begin
    a = AtomType("A"; wrap=false)
    b = AtomType("B"; wrap=false)
    # `wrap` is a keyword: the default wraps the SMARTS in `$(...)`, `wrap=false`
    # stores it verbatim.
    @test AtomType("A").smarts == "\$(A)"
    @test AtomType("A"; wrap=false).smarts == "A"
    @test_throws MethodError AtomType("A", false)
    # OR joins atom expressions with a comma; AND-NOT negates the second term.
    @test MG.smarts_or(a, b).smarts == "A,B"
    @test MG.smarts_andnot(a, b).smarts == "!B;A"
    # A raw SMARTS string on the right is wrapped as an AtomType first.
    @test MG.smarts_or(a, "B").smarts == MG.smarts_or(a, AtomType("B")).smarts
    # The deprecated +/- operators forward to the named combinators.
    @test (@test_deprecated a + b).smarts == MG.smarts_or(a, b).smarts
    @test (@test_deprecated a - b).smarts == MG.smarts_andnot(a, b).smarts
end

@testset "widened argument annotations" begin
    s = SubString("XA", 2)   # "A" as a SubString{String}
    # AtomType and the smarts combinators accept any AbstractString
    @test AtomType(s; wrap=false).smarts == "A"
    @test MG.smarts_or(AtomType("A"; wrap=false), s).smarts == "A,\$(A)"
    @test MG.smarts_andnot(AtomType("A"; wrap=false), s).smarts == "!\$(A);A"
    # FeatureDef accepts an AbstractString and an AbstractVector{<:Real}; its
    # fields are coerced to String / Vector{Float64} in the inner constructor.
    fd = MG.FeatureDef(s, :Carbon, view([1, 2, 3], 1:2))
    @test fd.smarts isa String && fd.smarts == "A"
    @test fd.weights isa Vector{Float64} && fd.weights == [1.0, 2.0]
    # feature_maps accepts AbstractVector arguments (views of fdefs / families)
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    carbondefs = [MG.FeatureDef("[#6]", :Carbon, [1.0])]
    @test haskey(feature_maps(mol, view(carbondefs, 1:1)), :Carbon)
    @test haskey(feature_maps(mol, FAMILY_DEFS, view([:Volume, :Donor], 1:1)), :Volume)
    # parse_feature_definitions accepts a SubString path
    defpath = joinpath(@__DIR__, "..", "assets", "const", "FeatureDefinitions.fdef")
    @test parse_feature_definitions(SubString(defpath)) isa MG.FamilyDef
end

@testset "malformed .fdef raises ArgumentError" begin
    # A negated atom type that references an undefined atom type is invalid input;
    # the parser reports it as an ArgumentError naming the offending type and line.
    mktemp() do path, io
        write(io, "AtomType !Undefined [#6]\n")
        close(io)
        @test_throws ArgumentError parse_feature_definitions(path)
        @test_throws "undefined atom type :Undefined" parse_feature_definitions(path)
    end
end

@testset ".fdef line continuation and bad lines" begin
    # Parse a .fdef written to a temporary file.
    parsefdef(content) = mktemp() do path, io
        write(io, content); close(io)
        return parse_feature_definitions(path)
    end
    # A trailing backslash continues a SMARTS onto the next line, which must be a
    # single whitespace-free token; the two pieces are concatenated.
    fdef = parsefdef("""
    AtomType Carbon [#6]\\
    [R]
    DefineFeature Cring [{Carbon}]\\
    [!H0]
      Family Rings
      Weights 1.0,2.0
    EndFeature
    """)
    @test MG.smarts(fdef.atomtypes[:Carbon]) == "[\$([#6][R])]"
    @test fdef.features[:Cring].family == :Rings
    @test fdef.features[:Cring].weights == [1.0, 2.0]
    # A continuation line carrying extra whitespace-separated tokens is rejected.
    @test_throws "continuation line must be a single token" parsefdef("AtomType Carbon [#6]\\\n[R] oops\n")
    # An unrecognized top-level keyword, and an unrecognized line inside a feature
    # block, are both reported as "bad input line for feature".
    @test_throws "bad input line for feature: BogusKeyword" parsefdef("BogusKeyword a b\n")
    @test_throws "bad input line for feature: BogusInner" parsefdef("DefineFeature F [#6]\n  BogusInner x\nEndFeature\n")
end

@testset "feature_maps merges shared families" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    # Two FeatureDefs assigned the same family are merged into one map entry whose
    # match list is the concatenation of the individual queries' matches.
    carbon = MG.FeatureDef("[#6]", :Heavy, [1.0])
    oxygen = MG.FeatureDef("[#8]", :Heavy, [1.0])
    merged = feature_maps(mol, [carbon, oxygen])
    @test collect(keys(merged)) == [:Heavy]
    ncarbon = length(feature_maps(mol, [carbon])[:Heavy])
    noxygen = length(feature_maps(mol, [oxygen])[:Heavy])
    @test length(merged[:Heavy]) == ncarbon + noxygen
end

@testset "compact show is newline-free" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    at = AtomType("A"; wrap=false)
    fd = MG.FeatureDef("[#6]", :Carbon, [1.0])
    pgmm = PharmacophoreGMM(mol)
    # the 2-arg show(io, x) is used to render container elements, so it must not
    # emit newlines — both alone and inside a vector.
    for x in (at, fd, pgmm)
        @test !occursin('\n', sprint(show, x))
        @test !occursin('\n', sprint(show, [x, x]))
    end
    # the pretty multi-line rendering lives on the MIME"text/plain" method.
    @test occursin('\n', sprint(show, MIME("text/plain"), pgmm))
    @test occursin('\n', sprint(show, MIME("text/plain"), fd))
end

@testset "MakieCore extension" begin
    # Loading MakieCore (at the top of this file) must activate the extension,
    # which supplies the pharmacophoredisplay method.
    ext = Base.get_extension(MolecularGaussians, :MolecularGaussiansMakieCoreExt)
    @test ext !== nothing
    @test !isempty(methods(MG.pharmacophoredisplay))
end

@testset "generic axes" begin
    # These functions annotate their array arguments as AbstractVector, promising
    # to work for any array. The result must depend only on the values, not on how
    # the collection is indexed: wrapping an input in an offset axis or a view must
    # leave the output unchanged. Their outputs (a feature-map Dict, a single
    # Gaussian) are not index-matched to any input, so only value-invariance is
    # asserted, not output axes.
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)

    # feature_maps(mol, fdefs): the vector of FeatureDefs is iterated, not indexed.
    fdefs = [MG.FeatureDef("[#6]", :Carbon, [1.0]), MG.FeatureDef("[#8]", :Oxygen, [1.0])]
    ref = feature_maps(mol, fdefs)
    @test feature_maps(mol, OffsetArray(fdefs, -3)) == ref
    @test feature_maps(mol, view(fdefs, :)) == ref

    # feature_maps(mol, familydef, families): the vector of family symbols likewise.
    families = [:Volume, :Donor]
    ref = feature_maps(mol, FAMILY_DEFS, families)
    @test feature_maps(mol, FAMILY_DEFS, OffsetArray(families, -3)) == ref
    @test feature_maps(mol, FAMILY_DEFS, view(families, :)) == ref

    # atoms_to_feature(mol, nodeset): the node-index collection is iterated to build
    # one Gaussian, for both the single-atom and multi-atom branches.
    nodes = collect(nodeset(mol))
    gausseq(a, b) = a.μ ≈ b.μ && a.σ ≈ b.σ && a.ϕ ≈ b.ϕ
    for ns in (nodes[1:1], nodes[1:4])
        ref = MG.atoms_to_feature(mol, ns)
        @test gausseq(MG.atoms_to_feature(mol, OffsetArray(ns, -3)), ref)
        @test gausseq(MG.atoms_to_feature(mol, view(ns, :)), ref)
    end
end

@testset "Aqua" begin
    # undocumented_names: docstring coverage is tracked separately.
    Aqua.test_all(MolecularGaussians; undocumented_names = false)
end

@testset "ExplicitImports" begin
    # `test_explicit_imports` also checks the MakieCore extension. The ignored
    # names have no public API in their defining packages: AbstractGMM,
    # AbstractStackedLabeledIsotropicGMM, ROCSAlignmentResult, centroid, distance,
    # local_align, tanimoto (GaussianMixtureAlignment) are the internals the core
    # builds on; @recipe, Theme, plot! (MakieCore) are the recipe interface and
    # Color, colortype (MolecularGraph) the color plumbing the extension builds
    # on; SimpleMolGraph (MolecularGraph) is the abstract supertype the core
    # dispatches on — the same coupling Aqua's piracy check exempts. The two
    # public-ness checks resolve "public" accurately only on Julia 1.11+, so they
    # are gated by version; the other five checks run everywhere.
    test_explicit_imports(MolecularGaussians;
        all_explicit_imports_are_public = VERSION >= v"1.11" ?
            (; ignore = (:AbstractGMM, :AbstractStackedLabeledIsotropicGMM, :ROCSAlignmentResult,
                         :centroid, :distance, :local_align, :tanimoto, Symbol("@recipe"),
                         :Theme, :plot!, :Color, :colortype, :SimpleMolGraph)) : false,
        all_qualified_accesses_are_public = VERSION >= v"1.11",
    )
end
