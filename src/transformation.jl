@static if (pkgversion(MolecularGraph) < v"0.22")
    getmass(a::SDFAtom) = a.mass
else
    getmass(a::SDFAtom) = a.isotope
end


"""
    transformed(tform, x)

Apply a coordinate transformation `tform` (any CoordinateTransformations
transform) to the atomic coordinates of `x`, returning a new `SDFAtom` or
`SDFMolGraph` whose coordinates are transformed and all other properties are
preserved.
"""
function transformed(tform, a::SDFAtom)
    newcoords = tform(a.coords)
    return SDFAtom(a.symbol, a.charge, a.multiplicity, getmass(a), Vector{eltype(newcoords)}(newcoords))
end

function transformed(tform, mol::SDFMolGraph)
    newmol = deepcopy(mol)
    for i in vertices(newmol)
        set_prop!(newmol, i, transformed(tform, props(mol, i)))
    end
    return newmol
end
