"""
    transformed(tform, x)

Apply a coordinate transformation `tform` (any CoordinateTransformations
transform) to the atomic coordinates of `x`, returning a new `SDFAtom` or
`SDFMolGraph` whose coordinates are transformed and all other properties are
preserved.
"""
function transformed(tform, a::SDFAtom)
    newcoords = tform(a.coords)
    return SDFAtom(a.symbol, a.charge, a.multiplicity, a.mass, newcoords)
end

function transformed(tform, mol::SDFMolGraph)
    newvprops = Dict(i => transformed(tform, a) for (i, a) in mol.vprops)
    return SDFMolGraph(mol.graph, newvprops, mol.eprops, mol.gprops, mol.state, mol.edge_rank)
end
