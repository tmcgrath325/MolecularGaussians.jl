import Base: eltype

# PharmacophoreGMM

struct PharmacophoreGMM{N,T<:Real,K} <: GaussianMixtureAlignment.AbstractLabeledIsotropicGMM{N,T,K}
    gaussians::Vector{IsotropicGaussian{N,T}}
    labels::Vector{K}
    axes::Vector{SVector{3,T}}           # bond axes (not necessarily normalized)
    origins::Vector{SVector{3,T}}        # bond origins (not necessarily centered on a gaussian)
    bondtofeatures::Vector{Vector{Int}} # gaussians that are rotated by a bond
    bondtobonds::Vector{Vector{Int}}     # other bonds that are rotated by a bond (further away from the molecule's "core")
end

eltype(::Type{PharmacophoreGMM{N,T,K}}) where {N,T,K} = GaussianMixtureAlignment.IsotropicGaussian{N,T}

"""
    model = PharmacophoreGMM(mol; bonds=nothing,
                                  combineatoms=false,
                                  rigid=false,
                                  σfun=vdwvolume_sigma,
                                  ϕfun=a->one(typeof(MolecularGraph.atom_radius(a))),
                                  featuremaps=Dict(:Volume => [[i] for i in heavy_atom_idxs(mol)]))
                            )

Creates a set Gaussian mixture models from a molecule or subgraph `mol`, with each model corresponding to a
particular type of molecular feature (e.g. hydrogen bond acceptors)

A mapping between pharmacophore feature types and vectors of atom indices, `featuremaps`, can be provided to
specify pharmacophore features. Otherwise, every atom will be assigned to individual `:Volume` features.
This mapping is typically generated as follows:
    ```
        fdefs = parse_feature_definitions()
        feats = [:Volume, :Hydrophobe, :Acceptor, :Donor, :PosIonizable, :NegIonizable] # this is not an exhaustive list
        feature_maps = feature_maps(mol, fdefs, feats)
    ```

The optional `combineatoms` keyword argument, if set to true, will generate single Gaussian features that combine grouped
atoms specified by `featuremaps` into single Gaussian features.

Optionally, functions `σfun` and `ϕfun` can be provided, which take an `SDFAtom` as input and return
variances `σ` and scaling coefficients `ϕ`, respectively.
By default, `σfun` is the `vdwvolume_sigma`, which corresponds to a value that ensure the volume of
the resulting Gaussian distribution matches the volume of the input atom as determined by its Van der Waals radius.
`ϕfun` defaults to scaling coefficients of 1.

The `bonds` optional keyword specifies bonds to be considered rotatable in the resulting model via a vector of
two-tuples containing atom indexes that define each bond. If not provided, all automatically detected rotatable bonds
will be used.

If `rigid` is provided and is true, no rotatable bonds will be given to the model, regardless of whether `bonds`
and `bondstoatoms` are provided.

"""
function PharmacophoreGMM(mol::SDFMolGraph;
                          bonds = nothing,
                          combineatoms = false,
                          rigid = false,
                          σfun = vdw_volume_sigma,
                          ϕfun = a -> one(typeof(MolecularGraph.atom_radius(a))),
                          featuremaps::Dict{K,Vector{Vector{Int}}} = Dict{Symbol,Vector{Vector{Int}}}(:Volume => [[i] for i in heavy_atom_idxs(mol)])) where K
    N = length(props(mol,1).coords)
    T = eltype(props(mol,1).coords)
    # prep for bond rotations
    sgs = rotatablesubgraphs(mol)
    if rigid
        bonds = Tuple{Int,Int}[]
        bondtoatoms = Vector{Int}[]
        bondtobonds = Vector{Int}[]
    elseif isnothing(bonds)
        bonds = [(sg.parentnodeidx, sg.childnodeidx) for sg in sgs]
    end
    bondtoatoms = [filter(x -> x ∉ bonds[i], sg.vlist) for (i,sg) in enumerate(sgs)]
    bondtobonds = [findall(b -> b[1] ∈ sg.vlist && b[2] ∈ sg.vlist, bonds) for sg in sgs]

    axes = [SVector{N,T}(props(mol,b[2]).coords .- props(mol,b[1]).coords) for b in bonds]
    origins = [SVector{N,T}(props(mol,b[1]).coords) for b in bonds]

    # add gaussians for each feature
    # TO DO: other/better options for keeping atoms separate (i.e. assigning multiple labels to a single atom)
    bondtofeatures = [Int[] for i=1:length(bonds)]
    gaussians = Vector{IsotropicGaussian{N,T}}()
    labels = Vector{Symbol}()
    for (feature, nodesets) in featuremaps
        for set in nodesets
            if combineatoms
                push!(gaussians, atoms_to_feature(mol, set; σfun=σfun, ϕfun=ϕfun))
                push!(labels, feature)
                for (i,atomsrotated) in enumerate(bondtoatoms)
                    len = length(filter(a -> a ∈ atomsrotated, set))
                    if (len >= length(set)/2)
                        push!(bondtofeatures[i], length(gaussians))
                    end
                end
            else
                for a in set
                    push!(gaussians, atoms_to_feature(mol, [a];  σfun=σfun, ϕfun=ϕfun))
                    for (i,atomsrotated) in enumerate(bondtoatoms)
                        if a ∈ atomsrotated
                            push!(bondtofeatures[i], length(gaussians))
                        end
                    end
                end
            end
        end
    end
    return PharmacophoreGMM{N,T,K}(gaussians, labels, axes, origins, bondtofeatures, bondtobonds)
end

function  Base.:*(R::AbstractMatrix{W}, x::PharmacophoreGMM{N,V,K}) where {N,V,K,W}
    numtype = promote_type(V, W)
    return PharmacophoreGMM{N,numtype,K}([R*g for g in x.gaussians], copy(x.labels), [R*ax for ax in x.axes], [R*o for o in x.origins], copy(x.bondtofeatures), copy(x.bondtobonds))
end

function  Base.:+(x::PharmacophoreGMM{N,V,K}, T::AbstractVector{W}) where {N,V,K,W}
    numtype = promote_type(V, W)
    return PharmacophoreGMM{N,numtype,K}([g+T for g in x.gaussians], copy(x.labels), [ax.+T for ax in x.axes], [o.+T for o in x.origins], copy(x.bondtofeatures), copy(x.bondtobonds))
end

Base.:-(x::PharmacophoreGMM, T::AbstractVector) = x + (-T)

feature_labels(x::PharmacophoreGMM) = unique(x.labels)

# descriptive display

Base.show(io::IO, pgmm::PharmacophoreGMM) = println(io,
    summary(pgmm),
    " with $(length(pgmm)) Gaussians with labels:\n",
    "$(feature_labels(pgmm))"
)
