import CoordinateTransformations 
import LinearAlgebra: I

## apply AffineMap to part of a graph

function (tform::AffineMap)(x::GraphMol{SDFileAtom, SDFileBond}, nodeset)
    return graphmol(x.edges, [i∈nodeset ? tform(n) : n for (i,n) in enumerate(nodeattrs(x))], edgeattrs(x))
end

## rotate coordinates about a particular axis

function angleaxis_rotate_coords(angle, axis, origin, args...)
    tform = Translation(origin) ∘ LinearMap(AngleAxis(angle, axis...)) ∘ Translation(-origin)
    return tform(args...)
end

angleaxis_rotate_nodes(nodes, angle, axis, origin, 
    ) = [angleaxis_rotate_coords(angle, axis, origin, n.coords) for n in nodes]

angleaxis_rotate_graph(graph::UndirectedGraph, nodeset, angle, axis, origin,
    ) = angleaxis_rotate_coords(angle, axis, origin, graph, nodeset)

angleaxis_rotate_graph(subgraph::SubgraphView, angle, axis, origin,
    ) = angleaxis_rotate_coords(angle, axis, origin, subgraph.graph, subgraph.nodes)


##  rotation about an edge of the graph

function rotated_components(graph::UndirectedGraph, edgeidx)
    # generate a subgraph after removing the specified edge, and obtain the nodes in each connected component
    edgeidxs = vcat(collect(1:edgeidx-1), collect(edgeidx+1:length(graph.edges)))
    nodesets = connectedcomponents(edgesubgraph(graph, edgeidxs))
    if length(nodesets) != 2
        throw(ArgumentError("Removing the specified edge does not generate two disjoint subgraphs."))
    end

    # identify which side of the edge will be rotated (the side with fewer nodes)
    sort!(nodesets; by=length)
    return nodesubgraph(graph, nodesets[1])
end

function rotate_edge(distalgraph::SubgraphView, proximalnode::SDFileAtom, distalnode::SDFileAtom, angle)
    axis = distalnode.coords .- proximalnode.coords
    origin = proximalnode.coords
    return angleaxis_rotate_graph(distalgraph, angle, axis, origin)
end

function rotate_edge(graph::UndirectedGraph, edgeidx, angle)
    distalgraph = rotated_components(graph, edgeidx)
    axisnodes = graph.edges[edgeidx]

    # specify the axis via the two nodes on either side of the removed edge
    distalnode   = nodeattr(graph, axisnodes[findfirst(x->x∈distalgraph.nodes, axisnodes)])
    proximalnode = nodeattr(graph, axisnodes[findfirst(x->x∉distalgraph.nodes, axisnodes)])

    # perform the rotation
    return rotate_edge(distalgraph, proximalnode, distalnode, angle)
end

function rotate_edges(graph::UndirectedGraph, edgeidxs, angles)
    rotgraph = graph
    for (i,edge) in enumerate(edgeidxs)
        rotgraph = rotate_edge(rotgraph, edge, angles[i])
    end
    return rotgraph
end
