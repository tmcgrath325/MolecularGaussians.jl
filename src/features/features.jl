"""
    feat = atoms_to_feature(mol::SDFMolGraph, nodeset; ϕfun=rocs_amplitude, σfun=vdw_volume_sigma)

Combine the atoms of `mol` indexed by `nodeset` into a single non-directional
`IsotropicGaussian` feature centered at their centroid.

`ϕfun(atom)` gives an atom's amplitude and `σfun(atom, ϕ)` its width. For a
single-atom `nodeset` the feature takes `ϕ = ϕfun(atom)` and `σ = σfun(atom, ϕ)`
directly. For several atoms, `ϕ` is the mean of `ϕfun` over them and `σ` is
derived from their van der Waals radii via the combined sphere volume (`σfun` is
not consulted in this case).
"""
function atoms_to_feature(mol::SDFMolGraph, nodeset; ϕfun = rocs_amplitude, σfun = vdw_volume_sigma)
    if length(nodeset)==1
        atom = props(mol, only(nodeset))
        μ = atom.coords
        σ = σfun(atom, ϕfun(atom))
        ϕ = ϕfun(atom)
    else
        atoms = [props(mol, node) for node in nodeset]
        coordmat = hcat([a.coords for a in atoms]...)
        μ = centroid(coordmat, fill(1/length(atoms), length(atoms)))
        ϕ = sum([ϕfun(a) for a in atoms])/length(atoms)
        σ = sphere_volume_sigma((sum(x -> x^3, [atom_radius(a) for a in atoms]))^(1/3), ϕ)
    end
    return IsotropicGaussian(μ, σ, ϕ)
end

# `public` is a soft keyword only on Julia 1.11+; guard so the module still
# loads on the 1.10 LTS, where the declaration is simply absent.
@static if VERSION >= v"1.11"
    eval(Meta.parse("public atoms_to_feature"))
end

function feature_maps(mol::SimpleMolGraph, fdefs::Vector{FeatureDef})
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

function feature_maps(mol::SimpleMolGraph, familydef::FamilyDef, families::Vector{Symbol}=collect(keys(familydef.families)))
    fdefs = reduce(vcat, [[familydef.features[name] for name in familydef.families[family]] for family in families])
    return feature_maps(mol, fdefs)
end