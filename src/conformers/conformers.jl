## filter out rotable bonds for terminal heavy atoms and hydrogens

function isterminalnode(mol, nodeidx, ignoreH=true)
    if ignoreH
        neighsymbols = [props(mol, n).symbol for n in neighbors(mol, nodeidx)]
        return count(x->x!=:H, neighsymbols) > 1 ? false : true
    else
        return length(neighbors(mol, nodeidx)) == 1 ? true : false
    end
end

function isterminalbond(mol, edge, ignoreH=true)
    term = (isterminalnode(mol, edge.src, ignoreH), isterminalnode(mol, edge.dst, ignoreH))
    return term[1] || term[2]
end

function rotatablebonds(mol::SDFMolGraph, ignoreH=true)
    rotatable = is_rotatable(mol)
    rotbonds = Vector{Graphs.SimpleEdge{Int}}()
    for (rot, e) in zip(rotatable, edges(mol))
        rot && !isterminalbond(mol, e, ignoreH) && push!(rotbonds, e)
    end
    # rank rotable bonds by how much effect they will have on shape
    sort!(rotbonds; rev=true, by=x->length(RotatableSubgraph(mol, x).vlist))
    return rotbonds
end

function rotablesubgraphs(mol::SDFMolGraph)
    return [RotatableSubgraph(mol, rb) for rb in rotatablebonds(mol)]
end



## generation of conformers for a particular bond-angle increment for all rotable bonds

function conformers(mol::SDFMolGraph; step=π/3, lower=-π, upper=π*(1-eps(Float64)), maxbonds=2)
    rotbonds = rotatablebonds(mol)
    if length(rotbonds) > maxbonds
        rotbonds = rotbonds[1:maxbonds]
    end
    anglevals = lower:step:upper
    confs = SDFMolGraph[]
    for i in CartesianIndices(tuple(fill(1:length(anglevals), length(rotbonds))...))
        angles = [anglevals[idx] for idx in Tuple(i)]
        push!(confs, rotate_edges(mol, rotbonds, angles))
    end
    return confs
end
