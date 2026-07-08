module MolecularGaussiansMakieCoreExt

using MolecularGaussians: MolecularGaussians, PharmacophoreGMM
import MakieCore: plot!
using MakieCore: @recipe, Theme
using MolecularGraph: stick!
using GaussianMixtureAlignment: gmmdisplay!, IsotropicGMM
using Colors: Colors

# Fallback palette cycled across GMM components that lack an entry in
# FEATURE_COLORS. CUD colors: https://jfly.uni-koeln.de/color/#assign
const DEFAULT_COLORS = [
    Colors.RGB(0   / 255, 114 / 255, 178 / 255), # blue
    Colors.RGB(230 / 255, 159 / 255, 0   / 255), # orange
    Colors.RGB(0   / 255, 158 / 255, 115 / 255), # green
    Colors.RGB(204 / 255, 121 / 255, 167 / 255), # reddish purple
    Colors.RGB(86  / 255, 180 / 255, 233 / 255), # sky blue
    Colors.RGB(213 / 255, 94  / 255, 0   / 255), # vermillion
    Colors.RGB(240 / 255, 228 / 255, 66  / 255), # yellow
]

const FEATURE_COLORS = Dict(
    :Donor         => Colors.RGB(255 / 255, 0   / 255, 255 / 255),  # magenta
    :Acceptor      => Colors.RGB(0   / 255, 255 / 255, 0   / 255),  # green
    :PosIonizable  => Colors.RGB(255 / 255, 0   / 255, 0   / 255),  # red
    :NegIonizable  => Colors.RGB(0   / 255, 0   / 255, 255 / 255),  # blue
    :Hydrophobe    => Colors.RGB(0   / 255, 255 / 255, 255 / 255),  # cyan
    :Ring          => Colors.RGB(255 / 255, 128 / 255, 255 / 255),  # orange
    :AromaticRing  => Colors.RGB(255 / 255, 64  / 255, 0   / 255),  # brown
    :Volume        => Colors.RGB(200 / 255, 200 / 255, 200 / 255),  # light grey
    :ExcludedVolume=> Colors.RGB(100 / 255, 100 / 255, 100 / 255),  # dark grey
)

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

function MolecularGaussians.pharmacophoredisplay(args...; kwargs...)
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
        allkeys = allkeys ∪ Set(mgmm.labels)
    end
    len = length(allkeys)
    for (i,k) in enumerate(allkeys)
        col = isnothing(color) ? (haskey(colors, k) ? colors[k] : palette[(i-1) % len + 1]) : color
        for mgmm in mgmms
            idxs = findall(isequal(k), mgmm.labels)
            isempty(idxs) || gmmdisplay!(md, IsotropicGMM(mgmm.gaussians[idxs]); display=disp, color=col, palette=palette)
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

end
