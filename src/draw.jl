using PlotlyJS

export drawMolGMM, drawMolGMMs, drawmol, drawPharmacophoreGMMs, plotdrawing

const atom_colors =     
    Dict(:C  => "#383838",   # dark grey
         :H  => "#b5b5b5",   # light grey
         :O  => "#d62728",   # red
         :N  => "#1f77b4",   # blue
         :S  => "#cbd123",   # yellow
         :Cl => "#2ca02c",   # green
         :P  => "#e29522")   # orange

const feature_colors = 
    Dict(:donor         => "#ff00ff",  # magenta
         :acceptor      => "#00ff00",  # green
         :cation        => "#ff0000",  # red
         :anion         => "#0000ff",  # blue
         :hydrophobe    => "#00ffff",  # cyan
         :rings         => "#ffa500",  # orange
         :aromaticrings => "#964b00")  # brown

function drawatoms(atoms::AbstractVector{SDFileAtom}; colordict=atom_colors, markersize=10)
    if isempty(atoms)
        return
    end
    color = colordict[atoms[1].symbol]  # all atoms assumed to be of the same element
    xs = [a.coords[1] for a in atoms]
    ys = [a.coords[2] for a in atoms]
    zs = [a.coords[3] for a in atoms]
    return scatter3d(;x=xs, y=ys, z=zs,
                     mode="markers", marker=attr(color=color, size=markersize, symbol="circle"),
                     showlegend=false, showscale=false, name="", hoverinfo="skip",
                     hovertemplate="")
end


function drawbonds(nodes::AbstractVector{SDFileAtom}, edges, edgeidxs, bonds::AbstractVector{SDFileBond}, planedir=SVector{3}([0.,0.,1.]); name = "bonds", linewidth=3, separation=0.2)
    xs, ys, zs = Float64[], Float64[], Float64[]
    for eidx in edgeidxs
        edge = edges[eidx]
        bondorder = bonds[eidx].order
        pos1, pos2 = nodes[edge[1]].coords, nodes[edge[2]].coords
        sepvec = cross(planedir, pos1.-pos2)
        sepvec = sepvec/norm(sepvec)
        if isodd(bondorder)
            d=0.
        else
            d=0.5
        end
        for i=1:Int(ceil(bondorder/2))
            dvec = (d+1-i)*separation*sepvec
            p1, p2 = pos1 .+ dvec, pos2 .+ dvec
            push!(xs, p1[1], p2[1], NaN)
            push!(ys, p1[2], p2[2], NaN)
            push!(zs, p1[3], p2[3], NaN)
            if (d+1-i)>0
                p1, p2 = pos1 .- dvec, pos2 .- dvec
                push!(xs, p1[1], p2[1], NaN)
                push!(ys, p1[2], p2[2], NaN)
                push!(zs, p1[3], p2[3], NaN)
            end
        end
    end

    return scatter3d(;x=xs, y=ys, z=zs,
                      mode="lines", line=attr(color="black", width=linewidth),
                      showlegend=!isnothing(name), name=isnothing(name) ? "" : name, hoverinfo="skip",
                      hovertemplate="")
end

drawbonds(mol::UndirectedGraph, planedir=SVector{3}([0.,0.,1.]); kwargs...
    ) = drawbonds(nodeattrs(mol), mol.edges, 1:length(mol.edges), edgeattrs(mol), planedir; kwargs...)

drawbonds(sg::SubgraphView, planedir=SVector{3}([0.,0.,1.]); kwargs...
    ) = drawbonds(nodeattrs(sg), sg.graph.edges, sg.edges, edgeattrs(sg), planedir; kwargs...)

function drawmol(mol::Union{GraphMol,UndirectedGraph,SubgraphView}; colordict=atom_colors, atoms=true, bonds=true, hydrogens=false, markersize=10, linewidth=1)
    traces = AbstractTrace[]
    if !hydrogens == true
        if typeof(mol) <: GraphMol
            mol = removehydrogens(mol)
        elseif typeof(mol) <: SubgraphView
            mol = nodesubgraph(removehydrogens(mol.graph), nodeset(removehydrogens(mol.graph)) ∩ nodeset(mol))
        end
    end

    formula = isa(mol,SubgraphView) ? molecularformula(mol.graph) : molecularformula(mol)
    elements = Symbol[]
    for idx in eachindex(formula)
        if isletter(formula[idx])
            # the last character should always be a number
            if idx < length(formula)
                if isletter(formula[idx+1])
                    push!(elements, Symbol(formula[idx:idx+1]))
                else
                    push!(elements, Symbol(formula[idx]))
                end
            else
                push!(elements, Symbol(formula[idx]))
            end
        end
    end

    coordmat = fill(NaN, 3, length(nodeset(mol)))
    # each atom type (CHONS, in most cases) gets its own trace to differentiate between colors, while keeping the number of traces low
    for el in elements
        elatoms = SDFileAtom[]
        for idx in nodeset(mol)
            if nodeattr(mol,idx).symbol == el
                coordmat[:,idx] = nodeattr(mol,idx).coords
                push!(elatoms, nodeattr(mol,idx))
            end
        end
        if atoms && !isempty(elatoms)
            push!(traces, drawatoms(elatoms; colordict=colordict, markersize=markersize))
        end
    end

    if bonds
        push!(traces, drawbonds(mol; linewidth=linewidth))
    end
    return traces
end

function drawMolGMM(molgmm::MolGMM; color=GaussianMixtureAlignment.default_colors[1], atoms=true, bonds=true, hydrogens=false, markersize=10, linewidth=1, kwargs...)
    sg = nodesubgraph(molgmm.graph, molgmm.nodes)
    traces = drawmol(sg; atoms=atoms, bonds=bonds, hydrogens=hydrogens, markersize=markersize, linewidth=linewidth)
    push!(traces, drawIsotropicGMM(molgmm; color=color, kwargs...)...)
    return traces
end

function drawMolGMMs(molgmms::AbstractVector{<:MolGMM}; colors=GaussianMixtureAlignment.default_colors, kwargs...)
    traces = AbstractTrace[]
    for (i,molgmm) in enumerate(molgmms)
        color = colors[mod(i-1, length(colors))+1]
        push!(traces, drawMolGMM(molgmm; color=color, kwargs...)...)
    end
    return traces
end

function drawPharmacophoreGMM(fgmm::PharmacophoreGMM; colordict=feature_colors, atoms=true, bonds=true, hydrogens=false, markersize=10, linewidth=1, kwargs...)
    sg = nodesubgraph(fgmm.graph, fgmm.nodes)
    traces = drawmol(sg; atoms=atoms, bonds=bonds, hydrogens=hydrogens, markersize=markersize, linewidth=linewidth)
    push!(traces, drawMultiGMM(fgmm; colordict=colordict, kwargs...)...)
    return traces
end

function drawPharmacophoreGMMs(fgmms::AbstractVector{<:PharmacophoreGMM}; colordict=feature_colors, kwargs...)
    traces = AbstractTrace[]
    for fgmm in fgmms
        push!(traces, drawPharmacophoreGMM(fgmm; colordict=colordict, kwargs...)...)
    end
    return traces
end
