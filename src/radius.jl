# Van der Waals radii in angstroms
const van_der_waals_radii = Dict{Int64, Float64}(MolecularGraph.ATOM_VANDERWAALS_RADII)

function vdwradii(mol::Union{UndirectedGraph,SubgraphView})
    return Dict(zip(nodeset(mol),[van_der_waals_radii[atomnumber(nodeattr(mol,idx).symbol)] for idx in nodeset(mol)]))
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

const const_volume_coeff = (2/(9π))^(1/3)

function vdwvolume_sigma(mol::Union{UndirectedGraph,SubgraphView}) 
    return Dict(zip(nodeset(mol),[const_volume_coeff * van_der_waals_radii[atomnumber(nodeattr(mol,idx).symbol)]^2 for idx in nodeset(mol)]))
end