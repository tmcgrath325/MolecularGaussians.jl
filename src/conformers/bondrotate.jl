import CoordinateTransformations 
import LinearAlgebra: I

## rotate about a particular axis, centered at a specified origin

function angleaxis_rotation(angle, axis, origin)
    return Translation(origin) ∘ LinearMap(RotMatrix(AngleAxis(angle, axis...))) ∘ Translation(-origin)
end

function angleaxis_rotate_coords(angle, axis, origin, args...)
    tform = angleaxis_rotation(angle, axis, origin)
    return tform(args...)
end

angleaxis_rotate_graph(graph::SDFMolGraph, angle, axis, origin,
    ) = angleaxis_rotate_coords(angle, axis, origin, graph)


##  rotation about an edge of the graph

struct RotatableSubgraph{T}
    edge::Graphs.SimpleEdge{Int}
    proximalidx::Int
    distalidx::Int
    vlist::Vector{Int}
    axis::SVector{3,T}
    origin::SVector{3,T}
end

function RotatableSubgraph(graph::SDFMolGraph, edge::Graphs.SimpleEdge)
    # generate a subgraph after removing the specified edge, and obtain the nodes in each connected component
    newgraph = deepcopy(graph)
    Graphs.rem_edge!(newgraph, edge)
    nodesets = connected_components(newgraph)
    if length(nodesets) != 2
        throw(ArgumentError("Removing the specified edge does not generate two disjoint subgraphs."))
    end

    # identify which side of the edge will be rotated (the side with fewer nodes)
    sort!(nodesets; by=length)
    reverseedge = edge.src ∈ nodesets[2]
    proximalidx = reverseedge ? edge.dst : edge.src
    distalidx = reverseedge ? edge.src : edge.dst
    
    proximalcoords = get_prop(graph, proximalidx, :coords)
    distalcoords = get_prop(graph, distalidx, :coords)

    T = eltype(distalcoords)
    axis = SVector{3,T}(distalcoords .- proximalcoords)
    origin = SVector{3,T}(distalcoords)
    subgraph, vlist = induced_subgraph(newgraph, nodesets[1])
    return RotatableSubgraph(edge, proximalidx, distalidx, vlist, axis, origin)
end

function rotate_edge!(graph::SDFMolGraph, edge, angle)
    rs = RotatableSubgraph(graph, edge)
    tform = angleaxis_rotation(angle, rs.axis, rs.origin)
    for i in rs.vlist
        graph.vprops[i] = tform(props(graph, i))
    end
    return graph
end
rotate_edge(graph, edge, angle) = rotate_edge!(deepcopy(graph), edge, angle)

function rotate_edges!(graph::SDFMolGraph, edgeidxs, angles)
    for (i,edge) in enumerate(edgeidxs)
        graph = rotate_edge!(graph, edge, angles[i])
    end
    return graph
end

rotate_edges(graph, edgeidxs, angles) = rotate_edges!(deepcopy(graph), edgeidxs, angles)
