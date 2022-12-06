# Van der Waals radii in angstroms
const van_der_waals_radii = Dict{Int64, Float64}(MolecularGraph.ATOM_VANDERWAALS_RADII)

function vdwradius(atom::SDFileAtom)
    return van_der_waals_radii[atomnumber(atom.symbol)]
end

function vdwradii(mol::Union{UndirectedGraph,SubgraphView})
    return Dict(zip(nodeset(mol), [vdwradius(nodeattr(mol,idx)) for idx in nodeset(mol)]))
end

"""
    vdwradius!(mol)

Add an attribute `:atomradius` to `mol`, a dictionary with the atom indices as keys and the Van der Waals 
radii as values.
"""
vdwradii!(mol::UndirectedGraph) = attributes(mol)[:vdwradii] = vdwradii(mol)

function vdwradii!(mol::SubgraphView)
    pmol = mol.graph
    pvdw = get!(attributes(pmol), :vdwradii, vdwradii!(pmol))
    return filter(pr -> pr.first ∈ nodeset(mol), pvdw)
end


const rocs_amplitude = 2.7
rocs_volume_amplitude(a) = rocs_amplitude
sphere_volume_sigma(r, ϕ) = (4/(3*ϕ*√π))^(1/3) * r

function vdw_volume_sigma(mol::Union{UndirectedGraph,SubgraphView}) 
    return Dict(zip(nodeset(mol),[const_volume_coeff * van_der_waals_radii[atomnumber(nodeattr(mol,idx).symbol)]^2 for idx in nodeset(mol)]))
end

vdw_volume_sigma(atom::SDFileAtom, ϕ = rocs_amplitude) = sphere_volume_sigma(vdwradius(atom), ϕ)