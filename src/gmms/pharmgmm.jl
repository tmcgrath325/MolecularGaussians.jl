mutable struct PharmacophoreGMM{N,T<:Real,K} <: AbstractIsotropicMultiGMM{N,T,K}
    gmms::Dict{K, FeatureGMM{N,T,K}}
    graph::Union{GraphMol{SDFileAtom, SDFileBond}, SubgraphView{GraphMol{SDFileAtom, SDFileBond}}}
    nodes::Set{Int}
    rotablesubgraphs::Vector{RotableSubgraph{N,T}}
    featuretypes::Vector{K}
    features::Dict{K,Vector{Set{Int}}}
end

eltype(::Type{PharmacophoreGMM{N,T,K}}) where {N,T,K} = FeatureGMM{N,T}

function combine(pgmm1::PharmacophoreGMM, pgmm2::PharmacophoreGMM)
    if dims(pgmm1) != dims(pgmm2)
        throw(ArgumentError("GMMs must have the same dimensionality"))
    end

    t = IsotropicGMM{dims(pgmm1),promote_type(numbertype(pgmm1), numbertype(pgmm2))}
    d = promote_type(typeof(pgmm1.gmms), typeof(pgmm2.gmms))
    gmms = d()
    xkeys, ykeys = keys(pgmm1.gmms), keys(pgmm2.gmms)
    for key in xkeys ∪ ykeys
        if key ∈ xkeys && key ∈ ykeys
            push!(gmms, Pair(key, convert(t, combine(pgmm1.gmms[key], pgmm2.gmms[key]))))
        elseif key ∈ xkeys
            push!(gmms, Pair(key, convert(t, pgmm1.gmms[key])))
        else
            push!(gmms, Pair(key, convert(t, pgmm2.gmms[key])))
        end
    end

    graph = deepcopy(pgmm1.graph)
    disjointunion!(graph, pgmm2.graph)

    rotsubgraphs = append!(pharmgmm1.rotablesubraphs, pharmgmm2.rotablesubraphs)
    sort!(rotsubgraphs, rev=true, by=x->length(x.nodes))

    firstlen = length(nodeattrs(pgmm1.graph))

    nodes = Set(pgmm1.nodes ∪ (pgmm2.nodes.+firstlen))

    featuretypes = pgmm1.featuretypes ∪ pgmm2.featuretypes

    features = Dict{eltype(xkeys ∪ ykeys), Vector{Set{Int}}}()
    for key in xkeys ∪ ykeys
        if key ∈ xkeys && key ∈ ykeys
            push!(features, key => vcat(pgmm1.features[key], [Set(st.+firstlen) for st in pgmm2.features[key]]))
        elseif key ∈ xkeys
            push!(features, key => pgmm1.features[key])
        else
            push!(features, key => [Set(st.+firstlen) for st in pgmm2.features[key]])
        end
    end

    return PharmacophoreGMM(gmms, graph, nodes, rotsubgraphs, featuretypes, features)
end

# rigid transformation

function transform!(pharmgmm::PharmacophoreGMM, tform::AffineMap)
    for (key, fgmm) in pharmgmm.gmms
        transform!(fgmm, tform)
    end
    pharmgmm.rotablesubgraphs = [tform(rsg) for rsg in pharmgmm.rotablesubgraphs]
    return pharmgmm
end

function (tform::AffineMap)(pharmgmm::PharmacophoreGMM)
    gmms = Dict([key => tform(fgmm) for (key, fgmm) in pharmgmm.gmms]...)
    return PharmacophoreGMM(gmms, pharmgmm.graph, pharmgmm.nodes, [tform(rsg) for rsg in pharmgmm.rotablesubgraphs], pharmgmm.featuretypes, pharmgmm.features)
end

# bond rotation

function bondrotate!(pharmgmm::PharmacophoreGMM, rotbondidx, angle)
    rotsubgraph = pharmgmm.rotablesubgraphs[rotbondidx]
    tform = angleaxis_rotation(angle, rotsubgraph.axis, rotsubgraph.origin)

    for (key,fgmm) in pharmgmm.gmms
        # bondrotate!(fgmm, rotbondidx, angle)
        push!(pharmgmm.gmms, key=>bondrotate(fgmm, rotsubgraph, tform))
    end
    pharmgmm.rotablesubgraphs[rotbondidx] = tform(rotsubgraph)
    return pharmgmm
end

function bondrotate(pharmgmm::PharmacophoreGMM{N,T,K}, rotbondidx, angle) where {N,T,K}
    n = N
    t = promote_type(T, typeof(angle))
    k = K

    rotsubgraph = pharmgmm.rotablesubgraphs[rotbondidx]
    tform = angleaxis_rotation(angle, rotsubgraph.axis, rotsubgraph.origin)

    gmms = Dict{k, FeatureGMM{n,t,k}}()
    for (key, fgmm) in pharmgmm.gmms
        push!(gmms, key=>bondrotate(fgmm, rotsubgraph, tform))
    end

    rotablesubgraphs = RotableSubgraph{n,t}[]
    for (i,rsg) in enumerate(pharmgmm.rotablesubgraphs)
        if i == rotbondidx
            push!(rotablesubgraphs, tform(rsg))
        else
            push!(rotablesubgraphs, rsg)
        end
    end

    return PharmacophoreGMM(gmms, pharmgmm.graph, pharmgmm.nodes, rotablesubgraphs, pharmgmm.featuretypes, pharmgmm.features)
end

function bondsrotate!(pharmgmm::PharmacophoreGMM, rotbondidxs, angles)
    for (i,ridx) in enumerate(rotbondidxs)
        bondrotate!(pharmgmm, ridx, angles[i])
    end
    return pharmgmm
end

function bondsrotate(pharmgmm::PharmacophoreGMM, rotbondidxs, angles)
    for (i,ridx) in enumerate(rotbondidxs)
        pharmgmm = bondrotate(pharmgmm, ridx, angles[i])
    end
    return pharmgmm
end

"""
    model = PharmacophoreGMM(mol, nodes=nodeset(mol); σfun=vdwvolume_sigma, ϕfun=ones, featurestypes=pubchem_features, features = pharmfeatures(mol), directional=true)

Creates a set Gaussian mixture models from a molecule or subgraph `mol`, with each model corresponding to a 
particular type of molecular feature (e.g. ring structures)

Optionally, functions `σfun` and `ϕfun` can be provided, which take `mol` as input and return dictionaries mapping
mapping node indicies to variances `σ` and scaling coefficients `ϕ`, respectively.

If `nodes` is provided, the Gaussian mixture models will be constructed only from atoms corresponding to the node
indexes of the molecule's graph.

`features` specifies the keys of allowed molecular features. Features with other keys will be ignored.

`directional` specifies whether or not geometric constraints for ring structures and hydrogen bond donors will be included.
"""
function PharmacophoreGMM(mol::UndirectedGraph,
                          nodes = nodeset(mol);
                          σfun = vdwvolume_sigma,
                          ϕfun = ones,
                          featuretypes = pubchem_features,
                          features = pharmfeatures(mol),
                          directional = true)
    dim = length(nodeattrs(mol)[1].coords)
    valtype = eltype(nodeattrs(mol)[1].coords)
    keytype = eltype(featuretypes)
    # rotable subgraphs
    rotsubgraphs = RotableSubgraph[]
    rotbonds = rotablebonds(mol)
    for rb in rotbonds
        push!(rotsubgraphs, RotableSubgraph(mol, rb))
    end
    # add a GMM for each type of feature
    gmms = Dict{Symbol,FeatureGMM{dim,valtype,keytype}}()
    for (key, nodesets) in features
        if key in featuretypes
            feats = FeatureGaussian{dim,valtype,keytype}[]
            for set in nodesets
                # check if all atoms in the feature are represented by allowed nodes
                if set ∩ nodes == set
                    # if it is a ring feature, add normal vectors
                    if (key == :rings || key == :aromaticrings) && directional
                        center = sum([nodeattr(mol,idx).coords for idx in set])/length(set)
                        coordmat = fill(zero(valtype), 3, length(set))
                        for (i,idx) in enumerate(set)
                            coordmat[:,i] = nodeattr(mol,idx).coords .- center
                        end
                        n = GaussianMixtureAlignment.planefit(coordmat)[1]
                        dirs = [n, -n]
                    # if it is a hydrogen bond donor, point outward from hydrogens
                    elseif key == :donor && directional
                        if length(set)==1 
                            dirs = Vector{valtype}[]   
                            idx = collect(set)[1]
                            for (edge,neigh) in neighbors(mol,idx)
                                if nodeattr(mol,neigh).symbol == :H
                                    push!(dirs, nodeattr(mol,neigh).coords - nodeattr(mol,idx).coords)
                                end                                
                            end
                        else
                            dirs = nothing
                        end
                    # if it is a hydrogen bond acceptor, point along the lone pair orbitals
                    elseif key == :acceptor && directional
                        # the acceptor should be a single atom
                        if length(set) != 1
                            dirs = nothing
                        else
                            idx = collect(set)[1]
                            neighs = collect(neighbors(mol,idx))
                            if hybridization(mol)[idx] == :sp3      # tetrahedral
                                if length(neighs) == 3
                                    # if there are 3 bonds
                                    bonddir1 = nodeattr(mol,neighs[1].second).coords .- nodeattr(mol,idx).coords
                                    bonddir2 = nodeattr(mol,neighs[2].second).coords .- nodeattr(mol,idx).coords
                                    bonddir3 = nodeattr(mol,neighs[3].second).coords .- nodeattr(mol,idx).coords
                                    dirs = Vector{valtype}[-normalize(bonddir1.+bonddir2.+bonddir3)]
                                elseif length(neighs) == 2
                                    # if there are 2 bonds
                                    bonddir1 = nodeattr(mol,neighs[1].second).coords .- nodeattr(mol,idx).coords
                                    bonddir2 = nodeattr(mol,neighs[2].second).coords .- nodeattr(mol,idx).coords
                                    axis1 = normalize(bonddir1.+bonddir2)
                                    axis2 = normalize(cross(bonddir1,bonddir2))
                                    R = AngleAxis(π,axis2...) * AngleAxis(π/2,axis1...)
                                    dirs = Vector{valtype}[R*bonddir1, R*bonddir2]
                                else
                                    dirs = nothing
                                end
                            elseif hybridization(mol)[idx] == :sp2  # trigonal planar
                                # if there are two bonds, use them to find the lone pair's direction
                                if length(neighs) != 2
                                    dirs = nothing
                                else
                                    dirs = Vector{valtype}[]
                                    bonddir1 = nodeattr(mol,neighs[1].second).coords .- nodeattr(mol,idx).coords
                                    bonddir2 = nodeattr(mol,neighs[2].second).coords .- nodeattr(mol,idx).coords
                                    axis = normalize(cross(bonddir1, bonddir2))
                                    push!(dirs, AngleAxis(4π/3, axis...)*bonddir1)
                                end
                            elseif hybridization(mol)[idx] == :sp   # linear
                                # direction points along the axis of the atom's single bond
                                if length(neighs) != 1
                                    dirs = nothing
                                else
                                    dirs = Vector{valtype}[nodeattr(nodeattr(mol,idx) .- mol,neighs[1].second.coords)]
                                end
                            end
                        end
                    # otherwise there is no geometric constraint
                    else
                        dirs = nothing
                    end
                    push!(feats, atoms_to_feature(mol, set, σfun, ϕfun, key, dirs))
                end
            end
            push!(gmms, Pair(key, FeatureGMM(feats, nodesets, key)))
        end
    end
    rotsubgraphs = rotablesubgraphs(mol)
    return PharmacophoreGMM(gmms, mol, nodes, rotsubgraphs, featuretypes, features)
end

PharmacophoreGMM(submol::SubgraphView; kwargs...) = PharmacophoreGMM(submol.graph, nodeset(submol); kwargs...)



Base.show(io::IO, pgmm::PharmacophoreGMM) = println(io,
    summary(pgmm),
    " from molecule with formula $(typeof(pgmm.graph)<:GraphMol ? molecularformula(pgmm.graph) : molecularformula(pgmm.graph.graph))",
    " with $(sum(length(gmm.second) for gmm in pgmm)) Gaussians in $(length(pgmm)) GMMs with labels:\n",
    "$([label for (label, gmm) in pgmm.gmms])"
)
