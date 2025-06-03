import MakieCore: plot!
using MakieCore: @recipe, Theme
using GaussianMixtureAlignment: gmmdisplay, DEFAULT_COLORS
using MolecularGraph: stick!
using Colors: coloralpha, RGB, N0f8

const FEATURE_COLORS = Dict(k => RGB{N0f8}((v ./ 255)...) for (k, v) in Dict(
    :Donor         => (255, 0,   255),  # magenta
    :Acceptor      => (0,   255, 0  ),  # green
    :PosIonizable  => (255, 0,   0  ),  # red
    :NegIonizable  => (0,   0,   255),  # blue
    :Hydrophobe    => (0,   255, 255),  # cyan
    :Ring          => (255, 128, 255),  # orange
    :AromaticRing  => (255, 64,  0  ),  # brown
    :Volume        => (200, 200, 200),  # light grey
    :ExcludedVolume=> (100, 100, 100),  # dark grey
))

@recipe(MolGMMDisplay, p) do scene
    Theme(
        display = :wire,
        color = nothing,
        palette = DEFAULT_COLORS,
        colors = FEATURE_COLORS,
        show_mol = true,
        mol_fun = stick!,
    )
end

function pharmacophoredisplay(args...; kwargs...)
    # set_theme!(backgroundcolor = :black)
    figaxplot = molgmmdisplay(args...; kwargs...)
    figaxplot.axis.show_axis[] = false
    return figaxplot
end

function plot!(md::MolGMMDisplay{<:NTuple{<:Any,<:PharmacophoreGMM{N,T,K}}}) where {N,T,K}
    mgmms = [md[i][] for i=1:length(md)]
    disp = md[:display][]
    color = md[:color][]
    colors = md[:colors][]
    palette = md[:palette][]
    allkeys = Set{K}()
    for mgmm in mgmms
        allkeys = allkeys ∪ keys(mgmm.gmms)
    end
    len = length(allkeys)
    for (i,k) in enumerate(allkeys)
        col = isnothing(color) ? (haskey(colors, k) ? colors[k] : palette[(i-1) % len + 1]) : color
        col  = isa(col, Color) ? coloralpha(col) : col
        for mgmm in mgmms
            haskey(mgmm.gmms, k) && gmmdisplay!(md, mgmm.gmms[k]; display=disp, color=col, palette=palette)
        end
    end
    if md[:show_mol][]
        molfun = md[:mol_fun][]
        for gmm in mgmms
            molfun(md, gmm.graph)
        end
    end
    return md
end