"""
    AtomIndexed(f)

Wrap a function so the feature-building machinery calls it with the molecule and the atom's
vertex index — `f(mol, i)` for an amplitude, `f(mol, i, ϕ)` for a width — instead of with the
atom's properties. Use for quantities that depend on the atom's graph context, such as
partial charges.
"""
struct AtomIndexed{F}
    f::F
end

# feature amplitude and width of one atom: plain functions see the atom's properties,
# `AtomIndexed` functions see the molecule and vertex index
_atom_value(f, mol, node) = f(props(mol, node))
_atom_value(f::AtomIndexed, mol, node) = f.f(mol, node)
_atom_value(f, mol, node, ϕ) = f(props(mol, node), ϕ)
_atom_value(f::AtomIndexed, mol, node, ϕ) = f.f(mol, node, ϕ)

"""
    feat = atoms_to_feature(mol::SDFMolGraph, nodeset; ϕfun=rocs_volume_amplitude, σfun=vdw_volume_sigma)

Combine the atoms of `mol` indexed by `nodeset` into a single non-directional
`IsotropicGaussian` feature centered at their centroid.

`ϕfun(atom)` gives an atom's amplitude and `σfun(atom, ϕ)` its width. For a
single-atom `nodeset` the feature takes `ϕ = ϕfun(atom)` and `σ = σfun(atom, ϕ)`
directly. For several atoms, `ϕ` is the mean of `ϕfun` over them and `σ` is
derived from their van der Waals radii via the combined sphere volume (`σfun` is
not consulted in this case).
"""
function atoms_to_feature(mol::SDFMolGraph, nodeset; ϕfun = rocs_volume_amplitude, σfun = vdw_volume_sigma)
    if length(nodeset)==1
        node = only(nodeset)
        atom = props(mol, node)
        μ = atom.coords
        ϕ = _atom_value(ϕfun, mol, node)
        σ = _atom_value(σfun, mol, node, ϕ)
    else
        atoms = [props(mol, node) for node in nodeset]
        coordmat = hcat([a.coords for a in atoms]...)
        μ = centroid(coordmat, fill(1/length(atoms), length(atoms)))
        ϕ = sum([_atom_value(ϕfun, mol, node) for node in nodeset])/length(nodeset)
        # Accumulate the combined volume at the coordinates' precision: `vdw_radius`
        # reads a Float32 table, and summing cubes there leaves the width dependent on
        # the order the atoms are listed in.
        T = eltype(μ)
        σ = sphere_volume_sigma(sum(a -> T(vdw_radius(a))^3, atoms)^(1/3), ϕ)
    end
    return IsotropicGaussian(μ, σ, ϕ)
end

# `public` is a soft keyword only on Julia 1.11+; guard so the module still
# loads on the 1.10 LTS, where the declaration is simply absent.
@static if VERSION >= v"1.11"
    eval(Meta.parse("public atoms_to_feature"))
end

"""
    feature_maps(mol, fdefs::AbstractVector{<:FeatureDef}) -> Dict{Symbol, Vector{Vector{Int}}}
    feature_maps(mol, familydef::FamilyDef, families=keys(familydef.families)) -> Dict{Symbol, Vector{Vector{Int}}}

Map each pharmacophore family to the sets of atom indices in `mol` that realize
it. For every `FeatureDef`, the atoms of each distinct substructure match of its
SMARTS pattern become one index set, grouped under the feature's family; each
such set is what [`atoms_to_feature`](@ref) later collapses into a single
Gaussian.

The second form draws its feature definitions from a [`FamilyDef`](@ref),
restricted to the requested `families` (all families by default).
"""
function feature_maps(mol::SimpleMolGraph, fdefs::AbstractVector{<:FeatureDef})
    featuremaps = Dict{Symbol,Vector{Vector{Int}}}()
    for fdef in fdefs
        query = smartstomol(smarts(fdef))
        iter = substruct_matches(mol, query)
        matchkeys = Set{Int}[]
        matches = Vector{Vector{Int}}()
        for it in iter
            itkeys = Set(keys(it))
            duplicate = false
            for m in matchkeys
                if itkeys == m
                    duplicate = true
                    break
                end
            end
            if !duplicate
                push!(matchkeys, itkeys)
                push!(matches, collect(itkeys))
            end
        end
        if haskey(featuremaps, fdef.family)
            featuremaps[fdef.family] = vcat(featuremaps[fdef.family], matches)
        else
            push!(featuremaps, fdef.family => matches)
        end
    end
    return featuremaps
end

function feature_maps(mol::SimpleMolGraph, familydef::FamilyDef, families::AbstractVector{Symbol}=collect(keys(familydef.families)))
    fdefs = [familydef.features[name] for family in families for name in familydef.families[family]]
    return feature_maps(mol, fdefs)
end