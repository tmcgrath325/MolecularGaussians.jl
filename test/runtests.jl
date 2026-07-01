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

using Graphs: induced_subgraph, edges, vertices, connected_components, rem_edge!
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
    # the type-based eltype agrees with the instance-based one (no stray graph
    # parameter leaking into the IsotropicGMM element type)
    @test eltype(typeof(pgmm)) == eltype(pgmm) == Pair{Symbol, MG.IsotropicGMM{3,Float64}}
end

@testset "atoms_to_feature defaults" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    remove_hydrogens!(mol)
    nodes = collect(nodeset(mol))
    # The documented defaults (ϕfun=rocs_volume_amplitude, σfun=vdw_volume_sigma)
    # build an IsotropicGaussian for both a single-atom and a multi-atom nodeset.
    @test MG.atoms_to_feature(mol, nodes[1:1]) isa MG.IsotropicGaussian
    @test MG.atoms_to_feature(mol, nodes[1:3]) isa MG.IsotropicGaussian
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
end

@testset "coordinate transforms" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "assets", "data", "E1050_3d.sdf"))
    orig1 = copy(props(mol, 1).coords)
    tform = AffineMap(RotationVec(0.3, -0.2, 0.5), SVector(1.0, 2.0, 3.0))
    # `transformed` maps every atom's coordinates through the tform, preserves the
    # coordinate container type, and leaves the source molecule unmutated.
    mol2 = MG.transformed(tform, mol)
    @test all(props(mol2, i).coords ≈ tform(props(mol, i).coords) for i in keys(mol.vprops))
    @test typeof(props(mol2, 1).coords) == typeof(props(mol, 1).coords)
    @test props(mol, 1).coords == orig1
    # Transforming a PharmacophoreGMM carries the transform into its molecular graph.
    pgmm = PharmacophoreGMM(mol)
    pgmm2 = MG.transform(pgmm, tform)
    @test all(props(pgmm2.graph, i).coords ≈ tform(props(mol, i).coords) for i in keys(mol.vprops))
    # `transform` accepts non-affine transforms too (rotation-only, translation-only).
    lin = LinearMap(RotationVec(0.3, -0.2, 0.5))
    tr = Translation(SVector(1.0, 2.0, 3.0))
    pgmm_lin = MG.transform(pgmm, lin)
    pgmm_tr = MG.transform(pgmm, tr)
    @test all(props(pgmm_lin.graph, i).coords ≈ lin(props(mol, i).coords) for i in keys(mol.vprops))
    @test all(props(pgmm_tr.graph, i).coords ≈ tr(props(mol, i).coords) for i in keys(mol.vprops))
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

@testset "Aqua" begin
    # undocumented_names: docstring coverage is tracked separately.
    Aqua.test_all(MolecularGaussians; undocumented_names = false)
end

@testset "ExplicitImports" begin
    # `test_explicit_imports` also checks the MakieCore extension. The ignored
    # names have no public API in their defining packages: AbstractIsotropicMultiGMM,
    # ROCSAlignmentResult, centroid, distance, local_align, tanimoto
    # (GaussianMixtureAlignment) are the internals the core builds on; @recipe,
    # Theme, plot! (MakieCore) are the recipe interface and Color, colortype
    # (MolecularGraph) the color plumbing the extension builds on — the same
    # coupling Aqua's piracy check exempts. The two public-ness checks resolve
    # "public" accurately only on Julia 1.11+, so they are gated by version; the
    # other five checks run everywhere.
    test_explicit_imports(MolecularGaussians;
        all_explicit_imports_are_public = VERSION >= v"1.11" ?
            (; ignore = (:AbstractIsotropicMultiGMM, :ROCSAlignmentResult, :centroid,
                         :distance, :local_align, :tanimoto, Symbol("@recipe"),
                         :Theme, :plot!, :Color, :colortype)) : false,
        all_qualified_accesses_are_public = VERSION >= v"1.11",
    )
end
