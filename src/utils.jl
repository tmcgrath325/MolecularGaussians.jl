nodeset(mol::SimpleMolGraph) = Set(vertices(mol))

heavy_atom_idxs(mol::SimpleMolGraph) = filter(x -> props(mol, x).symbol != :H, vertices(mol))