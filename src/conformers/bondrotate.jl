## rotate about a particular axis, centered at a specified origin

function angleaxis_rotation(angle, axis, origin)
    return Translation(origin) ∘ LinearMap(RotMatrix(AngleAxis(angle, axis...))) ∘ Translation(-origin)
end

function angleaxis_rotate_coords(angle, axis, origin, args...)
    tform = angleaxis_rotation(angle, axis, origin)
    return transformed(tform, args...)
end

angleaxis_rotate_graph(graph::SDFMolGraph, angle, axis, origin,
    ) = angleaxis_rotate_coords(angle, axis, origin, graph)


##  rotation about an edge of the graph

struct RotatableSubgraph{T}
    edge::Graphs.Edge{Int}
    proximalidx::Int
    distalidx::Int
    vlist::Vector{Int}
    axis::SVector{3,T}
    origin::SVector{3,T}
end

function RotatableSubgraph(graph::SDFMolGraph, edge::Graphs.Edge)
    # remove the edge and locate the components that hold each of its endpoints.
    # Only the endpoints' own component is split; any other components of a
    # disconnected graph are untouched and irrelevant to this bond's rotation.
    newgraph = deepcopy(graph)
    Graphs.rem_edge!(newgraph, edge)
    comps = connected_components(newgraph)
    srci = findfirst(c -> edge.src ∈ c, comps)
    dsti = findfirst(c -> edge.dst ∈ c, comps)
    if srci == dsti
        # the endpoints remain connected: the bond lies in a ring
        throw(ArgumentError("Removing the specified edge does not generate two disjoint subgraphs."))
    end

    # identify which side of the edge will be rotated (the side with fewer nodes)
    nodesets = sort([comps[srci], comps[dsti]]; by=length)
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
        set_prop!(graph, i, transformed(tform, props(graph, i)))
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


## rotation of a PharmacophoreGMM about one of its rotatable bonds

"""
    bondrotate(pgmm::PharmacophoreGMM, angle, bondidx::Int) -> PharmacophoreGMM
    bondrotate(pgmm::PharmacophoreGMM, angles, bondidxs::AbstractVector{Int}) -> PharmacophoreGMM

Rotate `pgmm` about its `bondidx`-th rotatable bond by `angle` radians, returning
a new `PharmacophoreGMM`. The components the bond moves (`pgmm.bondtogaussians[bondidx]`)
and the frames of the bonds downstream of it (`pgmm.bondtobonds[bondidx]`) are
rotated about the bond's axis; every other component and bond frame is left in place.

The stored `graph` is a rigid reference to the input conformer and is not updated
by a bond rotation, so after `bondrotate` it no longer matches the moved components.

The second form applies a sequence of rotations, `angles[k]` about `bondidxs[k]`,
in order; each rotation acts on the result of the previous one.
"""
function bondrotate(pgmm::PharmacophoreGMM{N,T,L,K,G}, angle, axis, origin, rotatedgaussians, rotatedbonds) where {N,T,L,K,G}
    tform = angleaxis_rotation(angle, axis, origin)
    rot = RotMatrix(AngleAxis(angle, axis...))
    # Promote to the angle's type so ForwardDiff duals (used to differentiate an alignment
    # objective through the rotation) survive rather than being truncated back to `T`.
    S = promote_type(T, typeof(angle))
    newgaussians = StackedLabeledGaussian{N,S,L,K}[i ∈ rotatedgaussians ? tform(g) : g for (i,g) in enumerate(pgmm.gaussians)]
    # axes are directions (rotate by the linear part only); origins are points
    newaxes = SVector{N,S}[i ∈ rotatedbonds ? SVector{N,S}(rot*a) : a for (i,a) in enumerate(pgmm.axes)]
    neworigins = SVector{N,S}[i ∈ rotatedbonds ? SVector{N,S}(tform(o)) : o for (i,o) in enumerate(pgmm.origins)]
    return PharmacophoreGMM{N,S,L,K,G}(newgaussians, pgmm.graph, pgmm.σfun, pgmm.ϕfun,
                                       pgmm.feature_maps, newaxes, neworigins, pgmm.bondtogaussians, pgmm.bondtobonds)
end

bondrotate(pgmm::PharmacophoreGMM, angle, bondidx::Int) =
    bondrotate(pgmm, angle, pgmm.axes[bondidx], pgmm.origins[bondidx], pgmm.bondtogaussians[bondidx], pgmm.bondtobonds[bondidx])

function bondrotate(pgmm::PharmacophoreGMM, angles, bondidxs::AbstractVector{Int})
    length(angles) == length(bondidxs) ||
        throw(DimensionMismatch("number of angles ($(length(angles))) must match number of bond indices ($(length(bondidxs)))"))
    newpgmm = pgmm
    for (angle, bondidx) in zip(angles, bondidxs)
        newpgmm = bondrotate(newpgmm, angle, bondidx)
    end
    return newpgmm
end
