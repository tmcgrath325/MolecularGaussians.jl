struct MolGMM{T<:Real,N}
    model::IsotropicGMM{T,N}
    graph::Union{UndirectedGraph,SubgraphView}
    nodes::Set{Int}
    σfun
    ϕfun
end
Base.eltype(::Type{MolGMM{T,N}}) where T where N where M = T
Base.length(gmm::MolGMM) = length(gmm.model.gaussians)
Base.size(gmm::MolGMM{T,N}) where T where N where M = (length(gmm), N)
Base.size(gmm::MolGMM, idx) = size(gmm)[idx]

"""
    model = MolGMM(mol, σfun=ones, ϕfun=ones, nodes=nodeset(mol))

Creates a Gaussian mixture model from a molecule or subgraph `mol`, which consists of spherical Gaussian distributions 
with means `μ` equal to atom coordinates. 

Optionally, functions `σfun` and `ϕfun` can be provided, which take `mol` as input and return dictionaries mapping
mapping node indicies to variances `σ` and scaling coefficients `ϕ`, respectively.

If `nodes` is provided, the Gaussian mixture model will be constructed only from atoms corresponding to the node
indexes of the molecule's graph.
"""
function MolGMM(mol::UndirectedGraph,
                σfun = ones, 
                ϕfun = ones,
                nodes=nodeset(mol))
    dim = length(nodeattrs(mol)[1].coords)
    valtype = eltype(nodeattrs(mol)[1].coords)
    σdict, ϕdict = Dict{Int, valtype}(σfun(mol)), Dict{Int, valtype}(ϕfun(mol))
    atoms = Vector{IsotropicGaussian{valtype,dim}}([IsotropicGaussian(SVector{dim}(nodeattr(mol, idx).coords), σdict[idx], ϕdict[idx]) for idx in nodes])
    # remove atoms with ϕ = 0
    filter!(atom->atom.ϕ≠0,atoms)
    return MolGMM(IsotropicGMM(atoms), mol, nodes, σfun, ϕfun)
end

MolGMM(submol::SubgraphView, σfun=ones, ϕfun=ones) = MolGMM(submol.graph, σfun, ϕfun, nodeset(submol))

struct FeatureMolGMM{T<:Real,N}
    model::MultiGMM{T,N}
    directional::Bool
    graph::Union{UndirectedGraph,SubgraphView}
    nodes::Set{Int}
    σfun
    ϕfun
end
Base.eltype(::Type{FeatureMolGMM{T,N}}) where T where N = T
Base.length(fgmm::FeatureMolGMM) = sum([length(fgmm.model.gmms[key]) for key in keys(fgmm.model.gmms)])
Base.size(fgmm::FeatureMolGMM{T,N}) where T where N = (sum([length(fgmm.model.gmms[key]) for key in keys(fgmm.model.gmms)]), N)
Base.size(fgmm::FeatureMolGMM, idx) = size(fgmm)[idx]

"""
    model = FeatureMolGMM(mol, σfun=vdwvolume_sigma, ϕfun=ones, nodes=nodeset(mol), features=pubchem_features)

Creates a set Gaussian mixture models from a molecule or subgraph `mol`, with each model corresponding to a 
particular type of molecular feature (e.g. ring structures)

Optionally, functions `σfun` and `ϕfun` can be provided, which take `mol` as input and return dictionaries mapping
mapping node indicies to variances `σ` and scaling coefficients `ϕ`, respectively.

If `nodes` is provided, the Gaussian mixture models will be constructed only from atoms corresponding to the node
indexes of the molecule's graph.

`features` specifies the keys of allowed molecular features. Features with other keys will be ignored.

`directional` specifies whether or not geometric constraints for ring structures and hydrogen bond donors will be included.
"""
function FeatureMolGMM(mol::UndirectedGraph,
                       σfun = vdwvolume_sigma, # vdwradii!,
                       ϕfun = ones,
                       nodes = nodeset(mol);
                       features = pubchem_features,
                       directional = true)
    dim = length(nodeattrs(mol)[1].coords)
    valtype = eltype(nodeattrs(mol)[1].coords)
    # add a GMM for each type of feature
    gmms = Dict{Symbol,IsotropicGMM{valtype,dim}}()
    for p in pharmfeatures!(mol)
        if p.first in features
            feats = IsotropicGaussian{valtype, dim}[]
            for set in p.second
                # check if all atoms in the feature are represented by allowed nodes
                if set ∩ nodes == set
                    # if it is a ring feature, add normal vectors
                    if (p.first == :rings || p.first == :aromaticrings) && directional
                        center = sum([nodeattr(mol,idx).coords for idx in set])/length(set)
                        coordmat = fill(zero(valtype), 3, length(set))
                        for (i,idx) in enumerate(set)
                            coordmat[:,i] = nodeattr(mol,idx).coords .- center
                        end
                        n = GOGMA.planefit(coordmat)[1]
                        dirs = [n, -n]
                    # if it is a hydrogen bond donor, point outward from hydrogens
                    elseif p.first == :donor && directional
                        dirs = Vector{valtype}[]   
                        for idx in set
                            for neigh in neighbors(mol,idx)
                                if nodeattr(mol,neigh.second).symbol == :H
                                    push!(dirs, nodeattr(mol,neigh.second).coords - nodeattr(mol,idx).coords)
                                end
                            end
                        end
                    # if it is a hydrogen bond acceptor, point along the lone pair orbitals
                    elseif p.first == :acceptor && directional
                        # the acceptor should be a single atom
                        if length(set) != 1
                            dirs = nothing
                        else
                            idx = collect(set)[1]
                            if hybridization(mol)[idx] == :sp3      # tetrahedral
                                # TODO
                                # need at least 2 bonds, then find the other directions
                                dirs = nothing
                            elseif hybridization(mol)[idx] == :sp2  # planar
                                # TODO
                                # need to get the plane using hybridization of the atom's neighbors (?)
                                dirs = nothing
                            end
                        end
                    # otherwise there is no geometric constraint
                    else
                        dirs = nothing
                    end
                    push!(feats, atoms_to_feature(mol, set, σfun, ϕfun, identity, dirs))
                end
            end
            push!(gmms, Pair(p.first, IsotropicGMM(feats)))
        end
    end
    gmms = gmms
    return FeatureMolGMM(MultiGMM(gmms), directional, mol, nodes, σfun, ϕfun)
end

FeatureMolGMM(submol::SubgraphView, σfun=ones, ϕfun=ones, features=pubchem_features, directional=true) = 
    FeatureMolGMM(submol.graph, σfun, ϕfun, nodeset(submol); features=features,directional=directional)

# descriptive display

Base.show(io::IO, molgmm::MolGMM) = println(io,
    "MolGMM from molecule with formula $(typeof(molgmm.graph)<:GraphMol ? molecularformula(molgmm.graph) : molecularformula(molgmm.graph.graph))",
    " with $(length(molgmm)) Gaussians."
)

Base.show(io::IO, fmolgmm::FeatureMolGMM) = println(io,
    "FeatureMolGMM from molecule with formula $(typeof(fmolgmm.graph)<:GraphMol ? molecularformula(fmolgmm.graph) : molecularformula(fmolgmm.graph.graph))",
    " with $(length(fmolgmm)) Gaussians in $(length(fmolgmm.model)) GMMs with labels:\n",
    "$([label for (label, gmm) in fmolgmm.model.gmms])"
)
