## filter out rotable bonds for terminal heavy atoms and hydrogens

function isterminalnode(mol, nodeidx, ignoreH=true)
    if ignoreH
        neighsymbols = [props(mol, n).symbol for n in neighbors(mol, nodeidx)]
        return count(x->x!=:H, neighsymbols) > 1 ? false : true
    else
        return length(neighbors(mol, nodeidx)) == 1 ? true : false
    end
end

function isterminalbond(mol, edge, ignoreH=true)
    term = (isterminalnode(mol, edge.src, ignoreH), isterminalnode(mol, edge.dst, ignoreH))
    return term[1] || term[2]
end

function rotatablebonds(mol::SDFMolGraph, ignoreH=true)
    rotatable = is_rotatable(mol)
    rotbonds = Vector{Graphs.Edge{Int}}()
    for (rot, e) in zip(rotatable, edges(mol))
        rot && !isterminalbond(mol, e, ignoreH) && push!(rotbonds, e)
    end
    # rank rotable bonds by how much effect they will have on shape
    sort!(rotbonds; rev=true, by=x->length(RotatableSubgraph(mol, x).vlist))
    return rotbonds
end

function rotablesubgraphs(mol::SDFMolGraph)
    return [RotatableSubgraph(mol, rb) for rb in rotatablebonds(mol)]
end



## generation of conformers for a particular bond-angle increment for all rotable bonds

function conformers(mol::SDFMolGraph; step=π/3, lower=-π, upper=π*(1-eps(Float64)), maxbonds=2)
    rotbonds = rotatablebonds(mol)
    if length(rotbonds) > maxbonds
        rotbonds = rotbonds[1:maxbonds]
    end
    anglevals = lower:step:upper
    confs = SDFMolGraph[]
    for i in CartesianIndices(tuple(fill(1:length(anglevals), length(rotbonds))...))
        angles = [anglevals[idx] for idx in Tuple(i)]
        push!(confs, rotate_edges(mol, rotbonds, angles))
    end
    return confs
end


function align_conformers(xconfs::AbstractVector{<:G}, yconfs::AbstractVector{<:G}; alignfun = local_align, kwargs...) where G<:PharmacophoreGMM
    bestx = xconfs[1]
    besty = yconfs[1]
    bestxidx = 1
    bestyidx = 1
    bestolap = Inf
    tform = identity
    for (i,yconf) in enumerate(yconfs)
        x, ctform, xidx, min = align_conformers(xconfs, yconf; alignfun, kwargs...)
        if min < bestolap
            bestx = x
            besty = yconf
            bestxidx = xidx
            bestyidx = i
            bestolap = min
            tform = ctform
        end
    end
    return bestx, besty, tform, bestxidx, bestyidx, bestolap
end

# Normalize each alignment backend's result to `(overlap, tform)`. Dispatch on the
# result type, not on the align function's identity, so wrapped or anonymous align
# functions work as long as they return one of the known result shapes.
_overlap_tform(res::ROCSAlignmentResult) = (res.minimum, res.tform)
_overlap_tform(res::Tuple) = (res[1], res[2])            # local_align
_overlap_tform(res) = (res.upperbound, res.tform_params) # gogma_align / tiv_gogma_align

function align_conformers(confs::AbstractVector{<:G}, template::L; alignfun = local_align, kwargs...) where {G<:PharmacophoreGMM, L<:AbstractIsotropicMultiGMM}
    bestconf = confs[1]
    bestidx = 1
    bestolap = Inf
    tform = identity
    for (i,conf) in enumerate(confs)
        min, restform = _overlap_tform(alignfun(conf, template; kwargs...))
        if min < bestolap
            bestconf = conf
            bestidx = i
            bestolap = min
            tform = restform
        end
    end
    return bestconf, tform, bestidx, bestolap
end