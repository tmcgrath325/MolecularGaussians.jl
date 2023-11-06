function Base.:+(a::SDFAtom, T::AbstractVector)
    newcoords = a.coords .+ T
    return SDFAtom(a.symbol, a.charge, a.multiplicity, a.mass, newcoords)
end

function Base.:*(R::AbstractMatrix, a::SDFAtom)
    newcoords = R * a.coords
    return SDFAtom(a.symbol, a.charge, a.multiplicity, a.mass, newcoords)
end

function Base.:+(mol::SDFMolGraph, T::AbstractVector)
    newvprops = deepcopy(mol.vprops)    
    for (i, a) in newvprops
        push!(newvprops, i => a + T)
    end
    return SDFMolGraph(mol.graph, newvprops, mol.eprops, mol.gprops, mol.state, mol.edge_rank)
end

function Base.:*(R::AbstractMatrix, mol::SDFMolGraph)
    newvprops = deepcopy(mol.vprops)    
    for (i, a) in newvprops
        push!(newvprops, i => R * a)
    end
    return SDFMolGraph(mol.graph, newvprops, mol.eprops, mol.gprops, mol.state, mol.edge_rank)
end