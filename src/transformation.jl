function (tform::AffineMap)(a::SDFileAtom)
    return SDFileAtom(a.symbol, a.charge, a.multiplicity, a.mass, tform(a.coords), a.stereo)
end

function (tform::AffineMap)(mol::GraphMol)
    return GraphMol(mol.neighbormap, mol.edges, [tform(a) for a in nodeattrs(mol)], mol.edgeattrs, mol.cache, mol.attributes)
end
