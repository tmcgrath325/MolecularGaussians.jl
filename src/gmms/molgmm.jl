
mutable struct MolGMM{N,T<:Real} <: AbstractIsotropicGMM{N,T}
    gaussians::Vector{AtomGaussian{N,T}}
    graph::Union{GraphMol{SDFileAtom, SDFileBond}, SubgraphView{GraphMol{SDFileAtom, SDFileBond}}}
    nodes::Set{Int}
    rotablesubgraphs::Vector{RotableSubgraph{N,T}}
end

eltype(::Type{MolGMM{N,T}}) where {N,T} = AtomGaussian{N,T}

convert(::Type{MolGMM{N,T}}, m::MolGMM) where {N,T} = MolGMM([convert(AtomGaussian{N,T},g) for g in m.gaussians], m.graph, m.nodes, [convert(RotableSubgraph{N,T}, rsg) for rsg in m.rotablesubgraphs])
convert(::Type{IsotropicGMM{N,T}}, m::MolGMM) where {N,T} = IsotropicGMM([convert(IsotropicGaussian{N,T},g) for g in m.gaussians])

function combine(molgmm1::MolGMM, molgmm2::MolGMM)
    if dims(molgmm1) != dims(molgmm2)
        throw(ArgumentError("GMMs must have the same dimensionality"))
    end
    t = promote_type(eltype(molgmm2), eltype(molgmm2))
    gaussians = t[vcat(molgmm1.gaussians, molgmm2.gaussians)...]

    graph = deepcopy(molgmm1.graph)
    disjointunion!(graph, molgmm2.graph)

    rotsubgraphs = append!(molgmm1.rotablesubraphs, molgmm2.rotablesubraphs)
    sort!(rotsubgraphs, rev=true, by=x->length(x.nodes))

    nodes = Set(molgmm1.nodes ∪ (molgmm2.nodes.+length(nodeattrs(molgmm1.graph))))

    return MolGMM(gaussians, graph, nodes, rotsubgraphs)
end


"""
    model = MolGMM(mol, nodes=nodeset(mol); σfun=ones, ϕfun=ones)

Creates a Gaussian mixture model from a molecule or subgraph `mol`, which consists of spherical Gaussian distributions 
with means `μ` equal to atom coordinates. 

Optionally, functions `σfun` and `ϕfun` can be provided, which take `mol` as input and return dictionaries mapping
mapping node indicies to variances `σ` and scaling coefficients `ϕ`, respectively.

If `nodes` is provided, the Gaussian mixture model will be constructed only from atoms corresponding to the node
indexes of the molecule's graph.
"""
function MolGMM(mol::UndirectedGraph,
                nodes=nodeset(mol);
                σfun = vdwvolume_sigma, 
                ϕfun = ones)
    N = length(nodeattrs(mol)[1].coords)
    T = eltype(nodeattrs(mol)[1].coords)
    σdict, ϕdict = Dict{Int, T}(σfun(mol)), Dict{Int, T}(ϕfun(mol))
    atoms = [AtomGaussian(nodeattr(mol, i), i, σdict[i], ϕdict[i]) for i in nodes]
    # remove atoms with ϕ = 0
    filter!(atom->atom.ϕ≠0,atoms)
    # prepare rotable subgraphs
    rotsubgraphs = rotablesubgraphs(mol)
    return MolGMM{N,T}(atoms, mol, nodes, rotsubgraphs)
end

MolGMM(submol::SubgraphView; kwargs...) = MolGMM(submol.graph, nodeset(submol); kwargs...)


# rigid transformation

function transform!(molgmm::MolGMM, tform::AffineMap)
    for ag in molgmm.gaussians
        transform!(ag, tform)
    end
    molgmm.rotablesubgraphs = [tform(rsg) for rsg in molgmm.rotablesubgraphs]
end

function (tform::AffineMap)(molgmm::MolGMM)
    return MolGMM([tform(g) for g in molgmm.gaussians], molgmm.graph, molgmm.nodes, [tform(rsg) for rsg in molgmm.rotablesubgraphs])
end


# bond rotations

function bondrotate!(molgmm::MolGMM, rotbondidx::Int, angle)
    rotsubgraph = molgmm.rotablesubgraphs[rotbondidx]
    tform = angleaxis_rotation(angle, rotsubgraph.axis, rotsubgraph.origin)

    for ag in molgmm.gaussians
        if ag.nodeidx ∈ rotsubgraph.nodes
            transform!(ag, tform)
        end
    end

    molgmm.rotablesubgraphs[rotbondidx] = tform(rotsubgraph)
    return molgmm
end

function bondrotate(molgmm::MolGMM{N,T}, rotbondidx::Int, angle) where {N,T}
    n = N
    t = promote_type(T, typeof(angle))
    rotsubgraph = molgmm.rotablesubgraphs[rotbondidx]
    tform = angleaxis_rotation(angle, rotsubgraph.axis, rotsubgraph.origin)

    gaussians = AtomGaussian{n,t}[]
    for ag in molgmm.gaussians
        if ag.nodeidx ∈ rotsubgraph.nodeidxs
            push!(gaussians, tform(ag))
        else
            push!(gaussians, ag)
        end
    end

    rotablesubgraphs = RotableSubgraph{n,t}[]
    for (i,rsg) in enumerate(molgmm.rotablesubgraphs)
        if i == rotbondidx
            push!(rotablesubgraphs, tform(rsg))
        else
            push!(rotablesubgraphs, rsg)
        end
    end

    return MolGMM(gaussians, molgmm.graph, molgmm.nodes, rotablesubgraphs)
end

function bondsrotate!(molgmm::MolGMM, rotbondidxs, angles)
    for (i,ridx) in enumerate(rotbondidxs)
        bondrotate!(molgmm, ridx, angles[i])
    end
    return molgmm
end

function bondsrotate(molgmm::MolGMM, rotbondidxs, angles)
    for (i,ridx) in enumerate(rotbondidxs)
        molgmm = bondrotate(molgmm, ridx, angles[i])
    end
    return molgmm
end

# descriptive display

Base.show(io::IO, molgmm::MolGMM) = println(io,
    "MolGMM from molecule with formula $(typeof(molgmm.graph)<:GraphMol ? molecularformula(molgmm.graph) : molecularformula(molgmm.graph.graph))",
    " with $(length(molgmm)) Gaussians."
)