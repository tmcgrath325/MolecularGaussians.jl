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
    parentnodeidx::Int
    childnodeidx::Int
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
    parentnodeidx = reverseedge ? edge.src : edge.dst
    childnodeidx = reverseedge ? edge.dst : edge.src
    
    parentnodecoords = get_prop(graph, parentnodeidx, :coords)
    childnodecoords = get_prop(graph, childnodeidx, :coords)

    T = eltype(childnodecoords)
    axis = SVector{3,T}(childnodecoords .- parentnodecoords)
    origin = SVector{3,T}(childnodecoords)
    subgraph, vlist = induced_subgraph(newgraph, nodesets[1])
    return RotatableSubgraph(edge, parentnodeidx, childnodeidx, vlist, axis, origin)
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


function bondrotate(pgmm::PharmacophoreGMM{N,T,K}, angle, axis, origin, rotatedfeatures, rotatedbonds) where {N,T,K}
    tform = angleaxis_rotation(angle, axis, origin)
    newgaussians = [i ∈ rotatedfeatures ? tform(g) : g for (i,g) in enumerate(pgmm.gaussians)]
    newaxes = [i ∈ rotatedbonds ? tform(a) : a for (i,a) in enumerate(pgmm.axes)]
    neworigins = [i ∈ rotatedbonds ? tform(o) : o for (i,o) in enumerate(pgmm.origins)]
    return PharmacophoreGMM{N,T,K}(newgaussians, newaxes, neworigins, copy(pgmm.bondtogaussians), copy(pgmm.bondtobonds))
end

bondrotate(pgmm::PharmacophoreGMM, angle, bondidx::Int) = bondrotate(pgmm, angle, pgmm.axes[bondidx], pgmm.origins[bondidx], pgmm.bondtogaussians[bondidx], pgmm.bondtobonds[bondidx])
function bondrotate(pgmm::PharmacophoreGMM, angles, bondidxs::AbstractVector{Int})
    newpgmm = bondrotate(pgmm, first(angles), first(bondidxs))
    for (angle, bondidx) in zip(angles[2:end],bondidxs[2:end])
        newpgmm = bondrotate(newpgmm, angle, bondidx)
    end
    return newpgmm
end
