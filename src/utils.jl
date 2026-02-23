using Graphs: vertices

nodeset(mol::SimpleMolGraph) = Set(vertices(mol))

is_terminal_bond(mol::SimpleMolGraph, edge::Graphs.SimpleEdge) = length(neighbors(mol, edge.src)) == 1 || length(neighbors(mol, edge.dst)) == 1

heavy_atom_idxs(mol::SimpleMolGraph) = filter(x -> props(mol, x).symbol != :H, vertices(mol))