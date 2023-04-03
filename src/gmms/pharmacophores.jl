const pubchem_features = [:acceptor, :anion, :cation, :donor, :hydrophobe, :rings, :aromaticrings]

"""
    feat = atoms_to_feature(atoms, dirs=nothing)

Combines Gaussian distributions specified by `atoms` to create a feature repressented by a 
single Gaussian. The feature is centered at the center of mass of the distributions, and the
width `σ` of the feature is the average of the distance of individual features from the center 
added to their widths.

Geometric constraints for the feature can be specified by `dirs`. The returned feature has no
direction by default.
"""
function atoms_to_feature(mol::UndirectedGraph, nodeset; ϕfun = rocs_amplitude, σfun = vdw_volume_sigma)
    if length(nodeset)==1
        atom = nodeattr(mol, only(nodeset))
        μ = atom.coords
        σ = σfun(atom, ϕfun(atom))
        ϕ = ϕfun(atom)
    else
        atoms = [nodeattr(mol, node) for node in nodeset]
        coordmat = hcat([a.coords for a in atoms]...)
        atomweights = [standardweight(a)[1] for a in atoms]
        μ = centroid(coordmat, atomweights ./ sum(atomweights))
        ϕ = sum([ϕfun(a) for a in atoms])/length(atoms)
        σ = sphere_volume_sigma((sum(x -> x^3, [vdwradius(a) for a in atoms]))^(1/3), ϕ)
    end

    return IsotropicGaussian(μ, σ, ϕ)
end
atoms_to_feature(submol::SubgraphView; kwargs...) = atoms_to_feature(submol.graph, nodeset(submol); kwargs...)

#TODO: fit single Gaussian to points sampled from sum of the features?
function combine_features(feats::AbstractVector{<:AbstractIsotropicGaussian}, weights=[f.ϕ for f in feats])
    ϕ = sum(weights)
    μ = centroid(feats, weights)
    σ = sum([f.σ * weights[i]/ϕ for (i,f) in enumerate(feats)])
    return IsotropicGaussian(μ, σ, ϕ)
end
