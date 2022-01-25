import CoordinateTransformations 
import LinearAlgebra: I

## apply AffineMap to part of a graph

function (tform::AffineMap)(x::GraphMol{SDFileAtom, SDFileBond}, nodeset)
    return graphmol(x.edges, [i∈nodeset ? tform(n) : n for (i,n) in enumerate(nodeattrs(x))], edgeattrs(x))
end

## rotate about a particular axis, centered at a specified origin

function angleaxis_rotation(angle, axis, origin)
    return Translation(origin) ∘ LinearMap(AngleAxis(angle, axis...)) ∘ Translation(-origin)
end

function angleaxis_rotate_coords(angle, axis, origin, args...)
    tform = angleaxis_rotation(angle, axis, origin)
    return tform(args...)
end

angleaxis_rotate_nodes(nodes, angle, axis, origin, 
    ) = [angleaxis_rotate_coords(angle, axis, origin, n.coords) for n in nodes]

angleaxis_rotate_graph(graph::UndirectedGraph, nodeset, angle, axis, origin,
    ) = angleaxis_rotate_coords(angle, axis, origin, graph, nodeset)

angleaxis_rotate_graph(subgraph::SubgraphView, angle, axis, origin,
    ) = angleaxis_rotate_coords(angle, axis, origin, subgraph.graph, subgraph.nodes)



##  rotation about an edge of the graph

struct RotableSubgraph{N,T}
    edgeidx::Int
    proximalidx::Int
    distalidx::Int
    nodeidxs::Set{Int}
    subgraph::SubgraphView{GraphMol{SDFileAtom, SDFileBond}}
    axis::SVector{N,T}
    origin::SVector{N,T}
end

function RotableSubgraph(edgeidx, proximalidx, distalidx, nodidxs, subgraph, axis, origin)
    t = promote_type(typeof(axis), typeof(origin))
    return RotableSubgraph(edgeidx, proximalidx, distalidx, nodidxs, subgraph, t(axis), t(origin))
end

function RotableSubgraph(graph::UndirectedGraph, edgeidx::Int)
    # generate a subgraph after removing the specified edge, and obtain the nodes in each connected component
    edgeidxs = vcat(collect(1:edgeidx-1), collect(edgeidx+1:length(graph.edges)))
    nodesets = connectedcomponents(edgesubgraph(graph, edgeidxs))
    if length(nodesets) != 2
        throw(ArgumentError("Removing the specified edge does not generate two disjoint subgraphs."))
    end

    # identify which side of the edge will be rotated (the side with fewer nodes)
    sort!(nodesets; by=length)
    edgenodeidxs = graph.edges[edgeidx]
    proximalidx = only(setdiff(edgenodeidxs, nodesets[1]))
    distalidx = only(setdiff(edgenodeidxs, proximalidx))
    
    proximalcoords = nodeattr(graph, proximalidx).coords
    distalcoords = nodeattr(graph, distalidx).coords

    n = length(distalcoords)
    t = eltype(distalcoords)
    axis = SVector{n,t}(distalcoords .- proximalcoords)
    origin = SVector{n,t}(distalcoords)
    return RotableSubgraph(edgeidx, proximalidx, distalidx, nodesets[1], nodesubgraph(graph, nodesets[1]), axis, origin)
end

convert(::Type{RotableSubgraph{N,T}}, rsg::RotableSubgraph) where {N,T} = RotableSubgraph(rsg.edgeidx, rsg.proximalidx, rsg.distalidx, rsg.nodeidxs, rsg.subgraph, SVector{N,T}(rsg.axis), SVector{N,T}(rsg.origin))

function (tform::AffineMap)(rsg::RotableSubgraph)
    return RotableSubgraph(rsg.edgeidx, rsg.proximalidx, rsg.distalidx, rsg.nodeidxs, rsg.subgraph, tform(rsg.axis), tform.linear * rsg.origin)
end


function rotate_edge(distalgraph::SubgraphView, angle, axis, origin)
    return angleaxis_rotate_graph(distalgraph, angle, axis, origin)
end
rotate_edge(rs::RotableSubgraph, angle) = rotate_edge(rs.subgraph, angle, rs.axis, rs.origin)

function rotate_edge!(graph::UndirectedGraph, edgeidx, angle)
    rs = RotableSubgraph(graph, edgeidx)
    nodes = nodeattrs(graph)
    rotated_subgraph = rotate_edge(rs, angle)
    for nodeidx in nodeset(rotated_subgraph)
        nodes[nodeidx] = nodeattr(rotated_subgraph, nodeidx)
    end
    return GraphMol(graph.neighbormap, graph.edges, nodes, graph.edgeattrs, graph.cache, graph.attributes)
end
rotate_edge(graph, edgeidx, angle) = rotate_edge!(deepcopy(graph), edgeidx, angle)

function rotate_edges!(graph::UndirectedGraph, edgeidxs, angles)
    for (i,edge) in enumerate(edgeidxs)
        graph = rotate_edge!(graph, edge, angles[i])
    end
    return graph
end

rotate_edges(graph, edgeidxs, angles) = rotate_edges!(deepcopy(graph), edgeidxs, angles)
