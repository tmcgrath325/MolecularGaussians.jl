using MolecularGraph: atomnumber
import MolecularGraph: atom_radius

MolecularGraph.atom_radius(a::SDFAtom) = MolecularGraph.ATOM_VANDERWAALS_RADII[atomnumber(a.symbol)]

const rocs_amplitude = 2.7

rocs_volume_amplitude(a) = rocs_amplitude
sphere_volume_sigma(r, ϕ) = (4/(3*ϕ*√π))^(1/3) * r

MolecularGraph.ATOM_VANDERWAALS_RADII

vdw_volume_sigma(mol::SimpleMolGraph, ϕ = rocs_amplitude) = [sphere_volume_sigma(r, ϕ) * r^2 for r in atom_radius(mol)]
vdw_volume_sigma(atom::SDFAtom, ϕ = rocs_amplitude) = sphere_volume_sigma(atom_radius(atom), ϕ)