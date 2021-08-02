add_attributes!(mol, attrs::Dict{Symbol}) = (merge!(mol.attributes, attrs); return mol)
add_attributes!(mol, allattrs::Dict{Int}) = add_attributes!(mol, allattrs[parse(Int, mol.attributes[:PUBCHEM_COMPOUND_CID])])

attributes(mol::UndirectedGraph) = mol.attributes
attributes(mol::UndirectedGraph, key::Symbol) = mol.attributes[key]

Base.zeros(mol::Union{UndirectedGraph,SubgraphView}) = Dict(zip(nodeset(mol), zeros(length(nodeset(mol)))))
Base.ones(mol::Union{UndirectedGraph,SubgraphView}) = Dict(zip(nodeset(mol), ones(length(nodeset(mol)))))

const_dict(val, mol::Union{UndirectedGraph,SubgraphView}) = Dict(zip(nodeset(mol), fill(val, length(nodeset(mol)))))