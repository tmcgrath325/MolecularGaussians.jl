"""
    vdw_radius(a::SDFAtom)

Van der Waals radius of atom `a`, looked up from
`MolecularGraph.ATOM_VANDERWAALS_RADII` by atomic number.
"""
vdw_radius(a::SDFAtom) = MolecularGraph.ATOM_VANDERWAALS_RADII[atom_number(a.symbol)]

const rocs_amplitude = 2.7

rocs_volume_amplitude(a) = rocs_amplitude
sphere_volume_sigma(r, ϕ) = (4/(3*ϕ*√π))^(1/3) * r

vdw_volume_sigma(atom::SDFAtom, ϕ = rocs_amplitude) = sphere_volume_sigma(vdw_radius(atom), ϕ)