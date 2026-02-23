using MolecularGraph: atom_number
import MolecularGraph: atom_radius

MolecularGraph.atom_radius(a::SDFAtom) = MolecularGraph.ATOM_VANDERWAALS_RADII[atom_number(a.symbol)]

const rocs_amplitude = 2.7  # Grant & Pickup, J. Phys. Chem. 1995
const rocs_lamda = 4 * π / 3 / rocs_amplitude
const rocs_k = π / rocs_lamda^(2/3)

rocs_volume_amplitude(a) = rocs_amplitude
sphere_volume_sigma(r, ϕ) = r / √2 * (4/(3*ϕ*√π))^(1/3)
sigma_to_radius(σ, ϕ) = √2 * σ / (4/(3*ϕ*√π))^(1/3)

MolecularGraph.ATOM_VANDERWAALS_RADII

vdw_volume_sigma(atom::SDFAtom, ϕ = rocs_amplitude) = sphere_volume_sigma(atom_radius(atom), ϕ)