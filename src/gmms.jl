import Base: eltype

# PharmacophoreGMM

struct PharmacophoreGMM{N,T<:Real,K} <: GaussianMixtureAlignment.AbstractLabeledIsotropicGMM{N,T,K}
    gaussians::Vector{LabeledIsotropicGaussian{N,T,K}}
    axes::Vector{SVector{3,T}}           # bond axes (not necessarily normalized)
    origins::Vector{SVector{3,T}}        # bond origins (not necessarily centered on a gaussian)
    bondtogaussians::Vector{Vector{Int}} # gaussians that are rotated by a bond
    bondtobonds::Vector{Vector{Int}}     # other bonds that are rotated by a bond (further away from the molecule's "core")
end

eltype(::Type{PharmacophoreGMM{N,T,K}}) where {N,T,K} = GaussianMixtureAlignment.LabeledIsotropicGaussian{N,T,K}

"""
    model = PharmacophoreGMM(mol, σfun=vdwvolume_sigma, ϕfun=ones, nodes=nodeset(mol), features=pubchem_features)

Creates a set Gaussian mixture models from a molecule or subgraph `mol`, with each model corresponding to a
particular type of molecular feature (e.g. ring structures)

Optionally, functions `σfun` and `ϕfun` can be provided, which take `mol` as input and return dictionaries
mapping node indicies to variances `σ` and scaling coefficients `ϕ`, respectively.

If `nodes` is provided, the Gaussian mixture models will be constructed only from atoms corresponding to the node
indexes of the molecule's graph.

`features` specifies the keys of allowed molecular features. Features with other keys will be ignored.

`directional` specifies whether or not geometric constraints for ring structures and hydrogen bond donors will be included.
"""
function PharmacophoreGMM(mol::SDFMolGraph,
                          bonds = nothing,
                          bondtoatoms = nothing,
                          bondtobonds = nothing;
                          combineatoms = false,
                          rigid = false,
                          σfun = vdw_volume_sigma,
                          ϕfun = a -> one(typeof(MolecularGraph.atom_radius(a))),
                          featuremaps::Dict{K,Vector{Vector{Int}}} = Dict{Symbol,Vector{Vector{Int}}}(:Volume => [[i] for i in heavy_atom_idxs(mol)])) where K
    N = length(props(mol,1).coords)
    T = eltype(props(mol,1).coords)
    # prep for bond rotations
    if isnothing(bonds) && isnothing(bondtoatoms) && isnothing(bondtobonds)
        if !rigid
            sgs = rotablesubgraphs(mol)
            bonds = [(sg.parentnodeidx, sg.childnodeidx) for sg in sgs]
            bondtoatoms = [filter(x -> x ∉ bonds[i], sg.vlist) for (i,sg) in enumerate(sgs)]
            bondtobonds = [findall(b -> b[1] ∈ sg.vlist && b[2] ∈ sg.vlist, bonds) for sg in sgs]
        else 
            bonds = Tuple{Int,Int}[]
            bondtoatoms = Vector{Int}[]
            bondtobonds = Vector{Int}[]
        end
    elseif rigid 
        bonds = Tuple{Int,Int}[]
        bondtoatoms = Vector{Int}[]
        bondtobonds = Vector{Int}[]
    end
    axes = [SVector{N,T}(props(mol,b[2]).coords .- props(mol,b[1]).coords) for b in bonds]
    origins = [SVector{N,T}(props(mol,b[1]).coords) for b in bonds]
    
    # add gaussians for each feature
    # TO DO: other/better options for keeping atoms separate (i.e. assigning multiple labels to a single atom)
    bondtofeatures = [Int[] for i=1:length(bonds)]
    gaussians = Vector{LabeledIsotropicGaussian{N,T,Symbol}}()
    for (feature, nodesets) in featuremaps
        for set in nodesets
            if combineatoms
                push!(gaussians, atoms_to_feature(mol, set, feature; σfun=σfun, ϕfun=ϕfun))
                for (i,atomsrotated) in enumerate(bondtoatoms)
                    len = length(filter(a -> a ∈ atomsrotated, set))
                    if (len >= length(set)/2)
                        push!(bondtofeatures[i], length(gaussians))
                    end
                end
            else
                for a in set
                    push!(gaussians, atoms_to_feature(mol, [a], feature;  σfun=σfun, ϕfun=ϕfun))
                    for (i,atomsrotated) in enumerate(bondtoatoms)
                        if a ∈ atomsrotated
                            push!(bondtofeatures[i], length(gaussians))
                        end
                    end
                end
            end
        end
    end
    return PharmacophoreGMM{N,T,K}(gaussians, axes, origins, bondtofeatures, bondtobonds)
end

function  Base.:*(R::AbstractMatrix{W}, x::PharmacophoreGMM{N,V,K}) where {N,V,K,W}
    numtype = promote_type(V, W)
    return PharmacophoreGMM{N,numtype,K}([R*g for g in x.gaussians], [R*ax for ax in x.axes], [R*o for o in x.origins], x.bondtogaussians, x.bondtobonds)
end

function  Base.:+(x::PharmacophoreGMM{N,V,K}, T::AbstractVector{W}) where {N,V,K,W}
    numtype = promote_type(V, W)
    return PharmacophoreGMM{N,numtype,K}([g+T for g in x.gaussians], [ax.+T for ax in x.axes], [o.+T for o in x.origins], x.bondtogaussians, x.bondtobonds)
end

Base.:-(x::PharmacophoreGMM, T::AbstractVector) = x + (-T)

feature_labels(x::PharmacophoreGMM) = unique([g.label for g in x.gaussians])

# descriptive display

Base.show(io::IO, pgmm::PharmacophoreGMM) = println(io,
    summary(pgmm),
    " with $(length(pgmm)) Gaussians with labels:\n",
    "$([label for label in feature_labels(pgmm)])"
)
