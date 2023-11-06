"""
    feat = atoms_to_feature(atoms, dirs=nothing)

Combines Gaussian distributions specified by `atoms` to create a feature repressented by a 
single Gaussian. The feature is centered at the center of mass of the distributions, and the
width `σ` of the feature is the average of the distance of individual features from the center 
added to their widths.

Geometric constraints for the feature can be specified by `dirs`. The returned feature has no
direction by default.
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

function feature_maps(mol::SimpleMolGraph, fdefs::Vector{FeatureDef})
    featuremaps = Dict{Symbol,Vector{Vector{Int}}}()
    for fdef in fdefs
        query = smartstomol(smarts(fdef))
        iter = substruct_matches(mol, query)
        matchkeys = Vector{Base.KeySet{Int64, Dict{Int64, Int64}}}()
        matches = Vector{Vector{Int}}()
        for it in iter
            itkeys = keys(it)
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