import Base: eltype

# PharmacophoreGMM

"""
    PharmacophoreGMM

A stacked labeled isotropic Gaussian mixture built from a molecule. Each component is a
`StackedLabeledGaussian` centered on one set of atoms and carrying `L` feature slots, with a
width, an amplitude, and a family label per slot. One point can therefore model several
physical interactions at once — steric bulk, electrostatics, hydrogen bonding — each sized
and weighted independently.

Feature sets built from the same atoms share a component. Components realizing fewer than
`L` features are padded with amplitude-zero slots, which contribute nothing.

Overlap between two `PharmacophoreGMM`s sums over the slot pairs whose labels interact (by
default, only equal labels; see `overlap` for the `interactions` keyword), computing the
distance between a pair of components once rather than once per pair of features.

The molecular `graph` is kept alongside the components so the molecule can still be drawn
and its formula reported. `σfun`, `ϕfun`, and `feature_maps` record how the features were
built.

The model also carries the geometry of its rotatable bonds so it can be flexed without the
graph (see [`bondrotate`](@ref)). `axes[b]` and `origins[b]` are the direction and a point on
the `b`-th bond's rotation axis; `bondtogaussians[b]` lists the components that bond moves
and `bondtobonds[b]` the bonds downstream of it. These four vectors have equal length (one
entry per rotatable bond).
"""
struct PharmacophoreGMM{N,T<:Real,L,K,G<:SimpleMolGraph} <: AbstractStackedLabeledIsotropicGMM{N,T,L,K}
    gaussians::Vector{StackedLabeledGaussian{N,T,L,K}}
    graph::G
    σfun
    ϕfun
    feature_maps::Dict{K, Vector{Vector{Int}}}
    axes::Vector{SVector{N,T}}
    origins::Vector{SVector{N,T}}
    bondtogaussians::Vector{Vector{Int}}
    bondtobonds::Vector{Vector{Int}}
    function PharmacophoreGMM{N,T,L,K,G}(gaussians, graph, σfun, ϕfun, feature_maps,
                                         axes, origins, bondtogaussians, bondtobonds) where {N,T,L,K,G}
        length(axes) == length(origins) == length(bondtogaussians) == length(bondtobonds) ||
            throw(DimensionMismatch("bond-geometry vectors (axes, origins, bondtogaussians, bondtobonds) must have equal length"))
        return new{N,T,L,K,G}(gaussians, graph, σfun, ϕfun, feature_maps, axes, origins, bondtogaussians, bondtobonds)
    end
end

function PharmacophoreGMM(gaussians::AbstractVector{StackedLabeledGaussian{N,T,L,K}},
                          graph::G, σfun, ϕfun, feature_maps,
                          axes, origins, bondtogaussians, bondtobonds) where {N,T,L,K,G<:SimpleMolGraph}
    return PharmacophoreGMM{N,T,L,K,G}(gaussians, graph, σfun, ϕfun, feature_maps,
                                       axes, origins, bondtogaussians, bondtobonds)
end

eltype(::Type{PharmacophoreGMM{N,T,L,K,G}}) where {N,T,L,K,G} = StackedLabeledGaussian{N,T,L,K}

"""
    feature_labels(pgmm::PharmacophoreGMM)

Return the distinct feature labels present in `pgmm`, in first-appearance order.
"""
feature_labels(pgmm::PharmacophoreGMM) = unique(l for g in pgmm.gaussians for l in g.labels)

"""
    family_function(f, family)

Select the width or amplitude function that applies to a pharmacophore `family`. A plain `f`
applies to every family; an `AbstractDict` maps each family to its own function and raises a
`KeyError` for a family it does not cover.
"""
family_function(f, family) = f
family_function(f::AbstractDict, family) = f[family]

"""
    PharmacophoreGMM(mol; combineatoms=true, rigid=false,
                          σfun=vdw_volume_sigma, ϕfun=a->one(typeof(vdw_radius(a))),
                          feature_maps=…)

Build a `PharmacophoreGMM` from a molecule or subgraph `mol`: one labeled feature per
molecular-feature set named in `feature_maps`, with features sharing a set of atoms stacked
onto a single component.

`feature_maps` maps each feature family (a `Symbol`) to a list of node-index sets. By
default it places one single-atom `:Volume` feature on every heavy atom of `mol`;
`feature_maps` derives pharmacophore families (donors, acceptors, rings, …) from feature
definitions. With `combineatoms=true` (the default) each set is collapsed into a single
feature by `atoms_to_feature`; with `combineatoms=false` each atom of a set gets its own
feature, all sharing the set's family label.

The model's rotatable bonds are detected automatically and their geometry recorded so it can
be flexed by [`bondrotate`](@ref). Pass `rigid=true` to record no rotatable bonds.

`ϕfun(atom)` gives an atom's amplitude and `σfun(atom, ϕ)` its width; see `atoms_to_feature`
for how they combine over multi-atom feature sets. Passing a `Dict` mapping families to
functions instead of a single function sizes and weights each family separately, which is
what makes stacking several features on one atom worthwhile.

# Example

```jldoctest
julia> mol = sdftomol(joinpath(pkgdir(MolecularGaussians), "assets", "data", "E1050_3d.sdf"));

julia> PharmacophoreGMM(mol)
PharmacophoreGMM{3, Float64, 1, Symbol, SDFMolGraph} from molecule with formula C18H24O8S2 with 28 stacked Gaussians labeled:
[:Volume]
```
"""
function PharmacophoreGMM(mol::SDFMolGraph;
                          combineatoms = true,
                          rigid = false,
                          σfun = vdw_volume_sigma,
                          ϕfun = a -> one(typeof(vdw_radius(a))),
                          feature_maps::Dict{K,Vector{Vector{Int}}} = Dict{Symbol,Vector{Vector{Int}}}(:Volume => [[i] for i in heavy_atom_idxs(mol)])) where K
    # SDFMolGraph coordinates are 3-D Float64, so an empty molecule (e.g. a
    # side chain with no atoms) still has a well-defined dimension and eltype.
    if isempty(vertices(mol))
        N, T = 3, Float64
    else
        N = length(props(mol,1).coords)
        T = eltype(props(mol,1).coords)
    end

    # rotatable-bond frames and, per bond, the atoms and downstream bonds it moves
    sgs = rigid ? RotatableSubgraph{T}[] : rotablesubgraphs(mol)
    axes = SVector{N,T}[SVector{N,T}(sg.axis) for sg in sgs]
    origins = SVector{N,T}[SVector{N,T}(sg.origin) for sg in sgs]
    # a bond's endpoints lie on its rotation axis and so do not move; exclude them
    bondtoatoms = [filter(a -> a != sg.edge.src && a != sg.edge.dst, sg.vlist) for sg in sgs]
    bondtobonds = [findall(s -> s.edge.src ∈ sg.vlist && s.edge.dst ∈ sg.vlist, sgs) for sg in sgs]

    # Features built from the same atoms have the same centroid, so they stack onto one
    # component keyed by that atom set.
    pointof = Dict{Vector{Int},Int}()
    atomsets = Vector{Int}[]
    μs = SVector{N,T}[]
    σs = Vector{T}[]
    ϕs = Vector{T}[]
    labelss = Vector{K}[]
    for (feature, nodesets) in feature_maps
        σf, ϕf = family_function(σfun, feature), family_function(ϕfun, feature)
        for set in nodesets
            for nodes in (combineatoms ? (collect(set),) : ([a] for a in set))
                # Sorting canonicalizes both the key and the centroid: summing coordinates
                # in the order a match happened to report them would put equal atom sets
                # at means differing in the last bit, and they would not stack.
                atoms = sort(nodes)
                feat = atoms_to_feature(mol, atoms; σfun=σf, ϕfun=ϕf)
                p = get!(pointof, atoms) do
                    push!(atomsets, atoms)
                    push!(μs, feat.μ)
                    push!(σs, T[]); push!(ϕs, T[]); push!(labelss, K[])
                    length(μs)
                end
                push!(σs[p], feat.σ); push!(ϕs[p], feat.ϕ); push!(labelss[p], feature)
            end
        end
    end

    bondtogaussians = [Int[] for _ in sgs]
    for (b, moved) in enumerate(bondtoatoms)
        for (p, atoms) in enumerate(atomsets)
            # a feature follows the bond when at least half of its atoms move
            count(a -> a ∈ moved, atoms) >= length(atoms)/2 && push!(bondtogaussians[b], p)
        end
    end

    gaussians = isempty(μs) ? StackedLabeledGaussian{N,T,1,K}[] :
        StackedLabeledIsotropicGMM(μs, σs, ϕs, labelss; padσ=one(T)).gaussians
    return PharmacophoreGMM(gaussians, mol, σfun, ϕfun, feature_maps,
                            axes, origins, bondtogaussians, bondtobonds)
end

# GaussianMixtureAlignment's generic `*`/`+` on AbstractStackedLabeledIsotropicGMM
# return a bare StackedLabeledIsotropicGMM, dropping the graph and feature metadata.
# These methods keep the PharmacophoreGMM type and carry those fields through.
# Bond `axes` are directions and `origins` are points: rotation moves both, but
# translation shifts only the points.
function  Base.:*(R::AbstractMatrix{W}, x::PharmacophoreGMM{N,V,L,K,G}) where {N,V,L,K,G,W}
    numtype = promote_type(V, W)
    gaussians = StackedLabeledGaussian{N,numtype,L,K}[R*g for g in x.gaussians]
    axes = SVector{N,numtype}[R*a for a in x.axes]
    origins = SVector{N,numtype}[R*o for o in x.origins]
    return PharmacophoreGMM{N,numtype,L,K,G}(gaussians, x.graph, x.σfun, x.ϕfun, x.feature_maps,
                                             axes, origins, x.bondtogaussians, x.bondtobonds)
end

function  Base.:+(x::PharmacophoreGMM{N,V,L,K,G}, T::AbstractVector{W}) where {N,V,L,K,G,W}
    numtype = promote_type(V, W)
    gaussians = StackedLabeledGaussian{N,numtype,L,K}[g+T for g in x.gaussians]
    axes = SVector{N,numtype}[a for a in x.axes]
    origins = SVector{N,numtype}[o .+ T for o in x.origins]
    return PharmacophoreGMM{N,numtype,L,K,G}(gaussians, x.graph, x.σfun, x.ϕfun, x.feature_maps,
                                             axes, origins, x.bondtogaussians, x.bondtobonds)
end

Base.:-(x::PharmacophoreGMM, T::AbstractVector) = x + (-T)

"""
    transform(pgmm::PharmacophoreGMM, tform) -> PharmacophoreGMM

Apply the transformation `tform` to every component of `pgmm` and to its underlying
molecular graph, returning a new `PharmacophoreGMM`. `tform` is called as `tform(g)` on each
`StackedLabeledGaussian` and must map a set of 3-D points rigidly (e.g. an `AffineMap` from
CoordinateTransformations). Bond `origins` are mapped as points and `axes` as directions (the
transform's linear part), recovered as `tform(o + a) - tform(o)`.
"""
function transform(pgmm::PharmacophoreGMM, tform)
    gaussians = [tform(g) for g in pgmm.gaussians]
    origins = [tform(o) for o in pgmm.origins]
    axes = [tform(o + a) - tform(o) for (o, a) in zip(pgmm.origins, pgmm.axes)]
    return PharmacophoreGMM(gaussians, transformed(tform, pgmm.graph), pgmm.σfun, pgmm.ϕfun,
                            pgmm.feature_maps, axes, origins, pgmm.bondtogaussians, pgmm.bondtobonds)
end

# `public` is a soft keyword only on Julia 1.11+; guard so the module still
# loads on the 1.10 LTS, where the declaration is simply absent.
@static if VERSION >= v"1.11"
    eval(Meta.parse("public transform"))
end

# descriptive display

Base.show(io::IO, pgmm::PharmacophoreGMM) = print(io,
    summary(pgmm),
    " with $(length(pgmm.gaussians)) stacked Gaussians labeled $(feature_labels(pgmm))"
)

Base.show(io::IO, ::MIME"text/plain", pgmm::PharmacophoreGMM) = print(io,
    summary(pgmm),
    " from molecule with formula $(molecular_formula(pgmm.graph))",
    " with $(length(pgmm.gaussians)) stacked Gaussians labeled:\n",
    "$(feature_labels(pgmm))"
)
