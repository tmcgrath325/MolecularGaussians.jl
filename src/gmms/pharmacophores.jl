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
function atoms_to_feature(atoms::AbstractVector{<:AbstractIsotropicGaussian}, dirs=nothing)
    if length(atoms) == 1
        return atoms[1]
    end
    μ = center_of_mass(atoms)
    # σ = mean([norm(μ-a.μ) + a.σ for (i,a) in enumerate(atoms)])
    σ = sum([norm(μ.-a.μ) for (i,a) in enumerate(atoms)])/length(atoms)
    ϕ = sum([a.ϕ for a in atoms])/length(atoms)
    if isnothing(dirs)
        dirs = typeof(μ)[]
    end
    return IsotropicGaussian(μ, σ, ϕ, dirs)
end
atoms_to_feature(gmm::AbstractSingleGMM, dirs=nothing) = atoms_to_feature(gmm.gaussians, dirs)
atoms_to_feature(mol::UndirectedGraph, nodeset, σfun, ϕfun, dirs=nothing) = atoms_to_feature(MolGMM(nodesubgraph(mol,nodeset), σfun, ϕfun), dirs)
atoms_to_feature(submol::SubgraphView, σfun, ϕfun, dirs=nothing) = atoms_to_feature(submol.graph, nodeset(submol), σfun, ϕfun, dirs)

#TODO: fit single Gaussian to points sampled from sum of the features?
function combine_features(feats::AbstractVector{<:AbstractIsotropicGaussian}, weights=[f.ϕ for f in feats])
    ϕ = sum(weights)
    μ = center_of_mass(feats, weights)
    σ = sum([f.σ * weights[i]/ϕ for (i,f) in enumerate(feats)])
    return IsotropicGaussian(μ, σ, ϕ)
end
