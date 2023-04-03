
import MakieCore: plot!
using MakieCore: @recipe, Theme, mesh!
using GeometryBasics: Point3f0, Sphere, Cylinder
using LinearAlgebra: cross, normalize
using MolecularGraph: nodeattr, nodeattrs, edgeattrs, neighbors, ATOMSYMBOLMAP, ATOM_VANDERWAALS_RADII, ATOM_COVALENT_RADII, UndirectedGraph

export spacefilling, spacefilling!, ballstick, ballstick!, stick, stick!, wire, wire!

struct Color 
    r::Int
    g::Int
    b::Int
end

const DEFAULT_BALL_WIDTH = Float32(0.3)
const DEFAULT_BOND_WIDTH = Float32(0.1)
const ATOM_COLORS = Dict(
    :default => Color(255,20,147),
    :H => Color(255,255,255),
    :B => Color(0,255,0),
    :C => Color(200,200,200),
    :N => Color(143,143,255),
    :O => Color(240,0,0),
    :F => Color(218,165,32),
    :Si => Color(218,165,32),
    :P => Color(255,165,0),
    :S => Color(255,200,50),
    :Cl => Color(0,255,0),
    :Br => Color(165,42,42),
    :I => Color(160,32,240),
    :Ba => Color(255,165,0),
    :Fe => Color(255,165,0),
    :Na => Color(0,0,255),
    :Mg => Color(34,139,34),
    :Zn => Color(165,42,42),
    :Cu => Color(165,42,42),
    :Ni => Color(165,42,42),
    :Ca => Color(128,128,144),
    :Mn => Color(128,128,144),
    :Al => Color(128,128,144),
    :Ti => Color(128,128,144),
    :Cr => Color(128,128,144),
    :Ag => Color(128,128,144),
    :Au => Color(218,165,32),
    :Li => Color(178,34,34),
    :He => Color(255,192,203)
)

colortype(c::Color) = Colors.RGB{Float32}(c.r/255, c.g/255, c.b/255)

function atom_radius(mol::UndirectedGraph, idx::Int; radius::Union{<:Real,<:Dict{Int,<:Any}} = ATOM_VANDERWAALS_RADII)
    isa(radius, Real) && return radius
    an = ATOMSYMBOLMAP[string(nodeattr(mol,idx).symbol)]
    r = radius[an]
    radius != ATOM_COVALENT_RADII && return Float32(r)
    if !isa(r, Real)
        # Carbon and a few metals have multiple values
        if an == 6
            nbrs = neighbors(mol, idx)
            key = nbrs == 3 ? "Csp3" :
                nbrs == 2 ? "Csp2" : "Csp"
            r = r[key]
        else
            # For metals, it's safest to choose the smallest radius
            r = minimum(values(r))
        end
    end
    return Float32(r)
end

"""
    spacefilling(mol::UndirectedGraph; radius="van der Waals", tform=identity)

Represent `mol` as a space-filling (Calotte) model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. (3D SDF files can be downloaded from sites such as PubChem.) The two supported options for `radius` are
`"van der Waals"` and `"covalent"`; the former are available only for main-group elements, and the latter are available for
all.

This function requires that you load one of the backends of the Makie/WGLMakie/CairoMakie family.
"""
function spacefilling(args...; radius="van der Waals", kwargs...)
    figaxplot = moldisplay(args...; radius=(radius == "covalent" ? ATOM_COVALENT_RADII : ATOM_VANDERWAALS_RADII), showbonds=false, kwargs...)
    fig, ax, md = figaxplot
    ax.show_axis[] = false
    return figaxplot
end
spacefilling!(args...; radius="van der Waals", kwargs...) = 
    moldisplay!(args...; radius=(radius == "covalent" ? ATOM_COVALENT_RADII : ATOM_VANDERWAALS_RADII), showbonds=false, kwargs...)

"""
    ballstick(mol::UndirectedGraph; radius=0.3, bondwidth=0.1)

Represent `mol` as a ball-and-stick model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. 3D SDF files can be downloaded from sites such as PubChem.

`radius` optionally specifies the radii of the balls, in Angstroms.
`bondwidth` optionally specifies the radii of the sticks, in Angstroms.

This function requires that you load one of the backends of the Makie/WGLMakie/CairoMakie family.
"""
function ballstick(args...; radius=DEFAULT_BALL_WIDTH, bondwidth=DEFAULT_BOND_WIDTH, kwargs...)
    figaxplot = moldisplay(args...; radius=Float32(radius), bondwidth=Float32(bondwidth), multiplebonds=true, kwargs...)
    fig, ax, md = figaxplot
    ax.show_axis[] = false
    return figaxplot
end
ballstick!(args...; radius=DEFAULT_BALL_WIDTH, bondwidth=DEFAULT_BOND_WIDTH, kwargs...) = moldisplay!(args...; radius=Float32(radius), bondwidth=Float32(bondwidth), multiplebonds=true, kwargs...)

"""
    stick(mol::UndirectedGraph; size=0.3)

Represent `mol` as a stick model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. 3D SDF files can be downloaded from sites such as PubChem.

`size` optionally specifies the width of the sticks, in Angstroms.

This function requires that you load one of the backends of the Makie/WGLMakie/CairoMakie family.
"""
function stick(args...; size=DEFAULT_BALL_WIDTH, kwargs...)
    figaxplot = moldisplay(args...; radius=Float32(size), bondwidth=Float32(size), multiplebonds=false, kwargs...)
    fig, ax, md = figaxplot
    ax.show_axis[] = false
    return figaxplot
end
stick!(args...; size=DEFAULT_BALL_WIDTH, kwargs...) = moldisplay!(args...; radius=Float32(size), bondwidth=Float32(size), multiplebonds=false, kwargs...)

"""
    wire(mol::UndirectedGraph; size=0.1)

Represent `mol` as a wire-frame model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. 3D SDF files can be downloaded from sites such as PubChem.

`size` optionally specifies the width of the bonds, in Angstroms.

This function requires that you load one of the backends of the Makie/WGLMakie/CairoMakie family.
"""
function wire(args...; size=DEFAULT_BOND_WIDTH, kwargs...) 
    figaxplot = moldisplay(args...; bondwidth=Float32(size), showatoms=false, multiplebonds=true, kwargs...)
    fig, ax, md = figaxplot
    ax.show_axis[] = false
    return figaxplot
end
wire!(args...; size=DEFAULT_BOND_WIDTH, kwargs...) = moldisplay!(args...; bondwidth=Float32(size), showatoms=false, multiplebonds=true, kwargs...)


@recipe(MolDisplay, mol) do scene
    Theme(
        radius = DEFAULT_BALL_WIDTH,
        bondwidth = DEFAULT_BOND_WIDTH,
        colors = ATOM_COLORS,
        multiplebonds = true,
        showbonds = true,
        showatoms = true,
    )
end

function plot!(md::MolDisplay{<:NTuple{<:Any,<:UndirectedGraph}})
    mols = [md[i][] for i=1:length(md)]
    colors = md[:colors][]
    for mol in mols
        if md[:showatoms][]
            for i=1:length(nodeattrs(mol))
                drawatom!(md, mol, i; radius=md[:radius][], colors=colors)
            end
        end
        if md[:showbonds][]
            for i=1:length(mol.edges)
                drawbond!(md, mol, i; bondwidth=md[:bondwidth][], multiplebonds=md[:multiplebonds][], colors=colors)
            end
        end
    end
    return md
end

function drawatom!(f, mol::UndirectedGraph, idx::Int; radius=DEFAULT_BALL_WIDTH, colors=ATOM_COLORS, kwargs...)
    r = atom_radius(mol, idx; radius=radius)
    a = nodeattr(mol,idx)
    p = Point3f0(a.coords)
    col = haskey(colors, a.symbol) ? colors[a.symbol] : colors[:default]
    mesh!(f, Sphere(p, r); color=colortype(col), kwargs...)
    return f
end

function drawbond!(f, mol::UndirectedGraph, idx::Int; bondwidth=DEFAULT_BOND_WIDTH, multiplebonds=false, colors=ATOM_COLORS, kwargs...)
    order = multiplebonds ? edgeattrs(mol)[idx].order : 1
    (atomidx1, atomidx2) = mol.edges[idx]
    a1, a2 = nodeattr(mol, atomidx1), nodeattr(mol, atomidx2)
    col1 = haskey(colors, a1.symbol) ? colors[a1.symbol] : colors[:default]
    col2 = haskey(colors, a2.symbol) ? colors[a2.symbol] : colors[:default]
    elementmatch = a1.symbol == a2.symbol
    pos1, pos2 = a1.coords, a2.coords   
    dir = normalize(cross([0,0,1], pos2 .- pos1))
    dists = (bondwidth * 2.5) .* collect(-0.5 * (order-1): 0.5 * (order-1))
    for dist in dists
        dvec = dist * dir
        p1, p2 = Point3f0((pos1 .+ dvec)...), Point3f0((pos2 .+ dvec)...)
        if elementmatch
            cyl = Cylinder(p1, p2, bondwidth)
            mesh!(f, cyl; color=colortype(col1), kwargs...)
        else
            midpoint = 0.5 .* (p1 .+ p2)
            pm = Point3f0(midpoint...)
            cyl1 = Cylinder(p1, pm, bondwidth)
            cyl2 = Cylinder(p2, pm, bondwidth)
            mesh!(f, cyl1; color=colortype(col1), kwargs...)
            mesh!(f, cyl2; color=colortype(col2), kwargs...)
        end
    end
    return f
end