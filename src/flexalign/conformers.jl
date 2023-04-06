## filter out rotable bonds for terminal heavy atoms and hydrogens

function isterminalnode(mol, nodeidx)
    neighsymbols = [nodeattr(mol, n).symbol for n in values(neighbors(mol, nodeidx))]
    return count(x->x!=:H, neighsymbols) > 1 ? false : true
end

function isterminalbond(mol, edgeidx)
    nodes = mol.edges[edgeidx]
    term = (isterminalnode(mol,nodes[1]), isterminalnode(mol, nodes[2]))
    return term[1] || term[2]
end

function rotablebonds(mol::GraphMol)
    rotbonds = findall(x->x==1, isrotatable(mol))
    filter!(x->!isterminalbond(mol,x), rotbonds)

    # rank rotable bonds by how much effect they will have on shape
    sort!(rotbonds; rev=true, by=x->length(RotableSubgraph(mol, x).nodeidxs))
    return rotbonds
end
rotablebonds(sg::SubgraphView) = rotablebonds(sg.graph)

function rotablesubgraphs(mol::GraphMol)
    return [RotableSubgraph(mol, rb) for rb in rotablebonds(mol)]
end



## generation of conformers for a particular bond-angle increment for all rotable bonds

function conformers(mol::UndirectedGraph; step=π/3, lower=-π, upper=π*(1-eps(Float64)), maxbonds=2)
    rotbonds = rotablebonds(mol)
    if length(rotbonds) > maxbonds
        rotbonds = rotbonds[1:maxbonds]
    end
    anglevals = lower:step:upper
    confs = GraphMol{SDFileAtom, SDFileBond}[]
    for i in CartesianIndices(tuple(fill(1:length(anglevals), length(rotbonds))...))
        angles = [anglevals[idx] for idx in Tuple(i)]
        push!(confs, rotate_edges(mol, rotbonds, angles))
    end
    return confs
end

function conformers(gmm::Union{MolGMM, PharmacophoreGMM}; step=π/3, lower=-π, upper=π*(1-eps(Float64)), maxbonds=2)
    rotsubgraphs = gmm.rotablesubgraphs
    if length(rotsubgraphs) > maxbonds
        rotsubgraphs = rotsubgraphs[1:maxbonds]
    end
    n = length(rotsubgraphs)
    anglevals = lower:step:upper
    confs = typeof(gmm)[]
    for i in CartesianIndices(tuple(fill(1:length(anglevals), n)...))
        angles = [anglevals[idx] for idx in Tuple(i)]
        push!(confs, bondsrotate(gmm, collect(1:n), angles))
    end
    return confs
end

