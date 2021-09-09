const default_colors =     
   ["#1f77b4",  # muted blue
    "#ff7f0e",  # safety orange
    "#2ca02c",  # cooked asparagus green
    "#d62728",  # brick red
    "#9467bd",  # muted purple
    "#8c564b",  # chestnut brown
    "#e377c2",  # raspberry yogurt pink
    "#7f7f7f",  # middle gray
    "#bcbd22",  # curry yellow-green
    "#17becf"]  # blue-teal

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

function plotdrawing(traces::AbstractVector{AbstractTrace}; size=20)
    layout = Layout(autosize=false, width=800, height=600,
                    margin=attr(l=0, r=0, b=0, t=65),
                    scene=attr(
                        aspectmode="cube",
                        xaxis=attr(visible=false, range=[-size/2, size/2]),
                        yaxis=attr(visible=false, range=[-size/2, size/2]),
                        zaxis=attr(visible=false, range=[-size/2, size/2]))
                    )
    plt = plot(traces, layout)
    return plt
end

function plotdrawing(traces::AbstractVector{<:AbstractVector{<:AbstractTrace}}; size=100)
    tracesvec = AbstractTrace[]
    for trs in traces
        push!(tracesvec, trs...)
    end
    layout = Layout(autosize=false, width=800, height=600,
                    margin=attr(l=0, r=0, b=0, t=65),
                    scene=attr(
                        aspectmode="cube",
                        xaxis=attr(visible=false, range=[-size/2, size/2]),
                        yaxis=attr(visible=false, range=[-size/2, size/2]),
                        zaxis=attr(visible=false, range=[-size/2, size/2]))
                    )
    plt = plot(tracesvec, layout)
    return plt
end

function drawatom(atom::SDFileAtom, tform=identity, colordict=atom_colors; markersize=10)
    pos = tform(atom.coords)
    color = colordict[atom.symbol]
    return scatter3d(;x=[pos[1]], y=[pos[2]], z=[pos[3]],
                     mode="markers", marker=attr(color=color, size=markersize, symbol="circle"),
                     showlegend=false, showscale=false, name="", hoverinfo="skip",
                     hovertemplate="")
end

function drawbond(atom1::SDFileAtom, atom2::SDFileAtom, tform=identity, bondorder=1, planedir=SVector{3}([0.,0.,1.]); linewidth=3, separation=0.2)
    pos1 = tform(atom1.coords)
    pos2 = tform(atom2.coords)
    traces = AbstractTrace[]
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
        push!(traces, scatter3d(;x=[p1[1],p2[1]], y=[p1[2],p2[2]], z=[p1[3],p2[3]],
                            mode="lines", line=attr(color="black", width=linewidth),
                            showlegend=false, name="", hoverinfo="skip",
                            hovertemplate=""))
        if (d+1-i)>0
            p1, p2 = pos1 .- dvec, pos2 .- dvec
            push!(traces, scatter3d(;x=[p1[1],p2[1]], y=[p1[2],p2[2]], z=[p1[3],p2[3]],
                                    mode="lines", line=attr(color="black", width=linewidth),
                                    showlegend=false, name="", hoverinfo="skip",
                                    hovertemplate=""))
        end
    end
    return traces
end

function drawmol(mol::Union{GraphMol,UndirectedGraph,SubgraphView}, tform=identity, colordict=atom_colors; atoms=true, bonds=true, hydrogens=false, markersize=10, linewidth=1)
    traces = AbstractTrace[]
    if !hydrogens == true
        if typeof(mol) <: GraphMol
            mol = removehydrogens(mol)
        elseif typeof(mol) <: SubgraphView
            mol = nodesubgraph(removehydrogens(mol.graph), nodeset(removehydrogens(mol.graph)) ∩ nodeset(mol))
        end
    end
    coordmat = fill(NaN, 3, length(nodeset(mol)))
    if atoms
        for (i,idx) in enumerate(nodeset(mol))
            coordmat[:,idx] = tform(nodeattr(mol,idx).coords)
            push!(traces, drawatom(nodeattr(mol,idx), tform, colordict; markersize=markersize))
        end
    end
    if bonds
        for bond in mol.edges
            if typeof(mol) <: SubgraphView
                bond = mol.graph.edges[bond]
            end
            bondorder = edgeattr(mol, findedgekey(mol, bond...)).order
            push!(traces, drawbond(nodeattr(mol,bond[1]), nodeattr(mol,bond[2]), tform, bondorder; linewidth=linewidth)...)
        end
    end
    return traces
end

function draw_wireframe_sphere(pos, r, w; npts=11, color=default_colors[1], name=nothing, opacity=1., sizecoef=1.)
    ϕs = range(0, stop=2π, length=npts)
    θs = range(-π/2, stop=π/2, length=npts)
    # vertical
    v_xs = [[r * sizecoef * cos(θ) * sin(ϕ) + pos[1] for θ in θs] for ϕ in ϕs[1:end-1]]
    v_ys = [[r * sizecoef * cos(θ) * cos(ϕ) + pos[2] for θ in θs] for ϕ in ϕs[1:end-1]]
    v_zs = [[r * sizecoef * sin(θ) + pos[3] for θ in θs] for ϕ in ϕs[1:end-1]]
    # horizontal
    h_xs = [[r * sizecoef * cos(θ) * sin(ϕ) + pos[1] for ϕ in ϕs] for θ in θs]
    h_ys = [[r * sizecoef * cos(θ) * cos(ϕ) + pos[2] for ϕ in ϕs] for θ in θs]
    h_zs = [[r * sizecoef * sin(θ) + pos[3] for ϕ in ϕs] for θ in θs] 

    hover = isnothing(name) ? "skip" : nothing

    xs, ys, zs = Float64[], Float64[], Float64[]
    for i=1:length(v_xs)
        append!(xs, v_xs[i])
        push!(xs, NaN)
        append!(ys, v_ys[i])
        push!(ys, NaN)
        append!(zs, v_zs[i])
        push!(zs, NaN)
    end
    for j=1:length(h_xs)
        append!(xs, h_xs[j])
        push!(xs, NaN)
        append!(ys, h_ys[j])
        push!(ys, NaN)
        append!(zs, h_zs[j])
        push!(zs, NaN)
    end

    return scatter3d(;x=xs, y=ys, z=zs, 
                      mode="lines",
                      line=attr(color=color),
                      opacity=opacity,
                      showlegend=false, showscale=false, name=name,
                      hovertemplate="μ = " * string(pos) * "<br>σ = " * string(r) * "<br>ϕ = " *string(w), 
                      hoverinfo=hover,
            )
end

function drawgaussian(gauss::IsotropicGaussian, tform=identity; sizecoef=1., opacity=1., color=default_colors[1], kwargs...)
    # gaussian centered at μ
    pos = tform(gauss.μ)
    r = gauss.σ * sizecoef
    gtrace = draw_wireframe_sphere(pos, gauss.σ, gauss.ϕ; color=color, kwargs...)
    if length(gauss.dirs) < 1
        return (gtrace,)
    end
    # cones to represent geometric constraints
    if typeof(tform) == typeof(identity)
        dirs = gauss.dirs
    else
        dirs = [tform.linear*dir for dir in gauss.dirs]
    end
    cntrs = [pos + 1.25*r*dir/norm(dir) for dir in dirs]
    dtrace = cone(;x=[c[1] for c in cntrs], y=[c[2] for c in cntrs], z=[c[3] for c in cntrs], 
                   u=[d[1] for d in dirs], v=[d[2] for d in dirs], w=[d[3] for d in dirs],
                   colorscale=[[0,color],[1,color]], opacity=opacity,
                   sizemode="absolute", sizeref=0.25,
                   showlegend=false, showscale=false, name="", hoverinfo="skip")
    return (gtrace, dtrace)
end

function drawIsotropicGMM(gmm::IsotropicGMM, tform=identity; kwargs...)
    # set opacities with weight values
    # weights = [gauss.ϕ for gauss in gmm.gaussians]
    # opacities = weights/maximum(weights) * 0.25

    # add a trace for each gaussian
    traces = AbstractTrace[]
    for i=1:length(gmm)
        push!(traces, drawgaussian(gmm.gaussians[i], tform; kwargs...)...)
    end
    return traces
end

function drawIsotropicGMMs(gmms::AbstractVector{<:IsotropicGMM},
                     tforms=fill(identity,length(gmms)); colors=default_colors, kwargs...)
    traces = AbstractTrace[]
    for (i,gmm) in enumerate(gmms)
        # add traces for each GMM
        color = colors[mod(i-1, length(colors))+1]
        push!(traces, drawIsotropicGMM(gmm, tforms[i]; color=color, kwargs...)...)
    end
    return traces
end

function drawMultiGMM(mgmm::MultiGMM, tform=identity; colordict=Dict{Symbol, String}(), colors=default_colors, kwargs...)
    # add traces from each GMM
    i = 1
    traces = AbstractTrace[]
    for key in keys(mgmm.gmms)
        # assign a color if the Dict doesn't include the key for a feature
        if key ∉ keys(colordict)
            push!(colordict, Pair(key, colors[mod(i-1, length(colors))+1]))
            i += i
        end
        push!(traces, drawIsotropicGMM(mgmm.gmms[key], tform; color=colordict[key], kwargs...)...)
    end
    return traces
end

function drawMultiGMMs(mgmms::AbstractVector{<:MultiGMM},
                            tforms=fill(identity,length(mgmms)); colordict=Dict{Symbol, String}(), colors=default_colors, kwargs...)
    # get all keys across the feature GMMs
    allkeys = Set{Symbol}()
    for fgmm in mgmms
        allkeys = allkeys ∪ keys(fgmm.gmms)
    end

    # assign a color to each feature type 
    i = 1
    for key in allkeys
        if key ∉ keys(colordict)
            push!(colordict, Pair(key, colors[mod(i-1, length(colors))+1]))
            i += i
        end
    end

    # add traces from each MultiGMM
    traces = AbstractTrace[]
    for (i,mgmm) in enumerate(mgmms)
        push!(traces, drawMultiGMM(mgmm, tforms[i]; colordict, colors=colors, kwargs...)...)
    end
    return traces
end

function drawMolGMM(molgmm::MolGMM, tform=identity; color=default_colors[1], kwargs...)
    sg = nodesubgraph(molgmm.graph, molgmm.nodes)
    traces = drawmol(sg, tform)
    push!(traces, drawIsotropicGMM(molgmm.model, tform; color=color, kwargs...)...)
    return traces
end

function drawMolGMMs(molgmms::AbstractVector{<:MolGMM}, tforms=fill(identity,length(molgmms)); colors=default_colors, kwargs...)
    traces = AbstractTrace[]
    for (i,molgmm) in enumerate(molgmms)
        color = colors[mod(i-1, length(colors))+1]
        push!(traces, drawMolGMM(molgmm, tforms[i]; color=color, kwargs...)...)
    end
    return traces
end

function drawFeatureMolGMM(fgmm::FeatureMolGMM, tform=identity; colordict=feature_colors, kwargs...)
    sg = nodesubgraph(fgmm.graph, fgmm.nodes)
    traces = drawmol(sg, tform)
    push!(traces, drawMultiGMM(fgmm.model, tform; colordict=colordict, kwargs...)...)
    return traces
end

function drawFeatureMolGMMs(fgmms::AbstractVector{<:FeatureMolGMM}, tforms=fill(identity,length(fgmms));
                            colordict=feature_colors, kwargs...)
    traces = AbstractTrace[]
    for (i,fgmm) in enumerate(fgmms)
        push!(traces, drawFeatureMolGMM(fgmm, tforms[i]; colordict=colordict, kwargs...)...)
    end
    return traces
end
