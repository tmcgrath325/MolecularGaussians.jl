function (tform::AffineMap)(a::SDFileAtom)
    return SDFileAtom(a.symbol, a.charge, a.multiplicity, a.mass, tform(a.coords), a.stereo)
end

function (tform::AffineMap)(mol::GraphMol)
    return GraphMol(mol.neighbormap, mol.edges, [tform(a) for a in nodeattrs(mol)], mol.edgeattrs, mol.cache, mol.attributes)
end

# function (tform::AffineMap)(x::MolGMM)
#     T = typeof(x)
#     return T.name.wrapper([tform(g) for g in x.gaussians], tform(x.graph), x.nodes, x.σfun, x.ϕfun)
# end

# function (tform::AffineMap)(x::PharmacophoreGMM)
#     T = typeof(x)
#     tformgmms = [tform(x.gmms[key]) for key in keys(x.gmms)]
#     gmmdict = Dict{eltype(keys(x.gmms)),eltype(tformgmms)}()
#     for (i,key) in enumerate(keys(x.gmms))
#         push!(gmmdict, key=>tformgmms[i])
#     end
#     return T.name.wrapper(gmmdict, x.directional, tform(x.graph), x.nodes, x.σfun, x.ϕfun)
# end
