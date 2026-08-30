## Gasteiger–Marsili partial atomic charges by iterative partial equalization of orbital
## electronegativity (PEOE): Gasteiger & Marsili, Tetrahedron 36 (1980) 3219.
##
## Each atom's electronegativity is a quadratic in its charge, χ(q) = a + b·q + c·q², with
## (a, b, c) depending on element and hybridization. Every iteration moves charge across each
## bond in proportion to the electronegativity difference, scaled by the cation
## electronegativity χ⁺ = a + b + c of the donating (less electronegative) atom and damped by
## 2⁻ᵏ at iteration k, so the transfers converge geometrically. Parameter values follow the
## RDKit implementation, the de facto reference for these charges.

# (a, b, c) per (element, mode). The mode is the hybridization, except: hydrogen and the
# halogens have a single parameter set, and a sulfur whose hybridization is unperceived
# (:none) falls back to oxygen-count modes (:so, :so2) — matching the RDKit reference
# implementation, which uses the oxidized-sulfur parameters only in that case.
const GASTEIGER_PARAMS = Dict{Tuple{Symbol,Symbol},NTuple{3,Float64}}(
    (:H, :any)  => (7.17, 6.24, -0.56),
    (:C, :sp3)  => (7.98, 9.18, 1.88),
    (:C, :sp2)  => (8.79, 9.32, 1.51),
    (:C, :sp)   => (10.39, 9.45, 0.73),
    (:N, :sp3)  => (11.54, 10.82, 1.36),
    (:N, :sp2)  => (12.87, 11.15, 0.85),
    (:N, :sp)   => (15.68, 11.70, -0.27),
    (:O, :sp3)  => (14.18, 12.92, 1.39),
    (:O, :sp2)  => (17.07, 13.79, 0.47),
    (:F, :any)  => (14.66, 13.85, 2.31),
    (:Cl, :any) => (11.00, 9.69, 1.35),
    (:Br, :any) => (10.08, 8.47, 1.16),
    (:I, :any)  => (9.90, 7.96, 0.96),
    (:S, :sp3)  => (10.14, 9.13, 1.38),
    (:S, :sp2)  => (10.88, 9.49, 1.33),
    (:S, :so)   => (10.14, 9.13, 1.38),
    (:S, :so2)  => (12.00, 10.81, 1.20),
    (:P, :sp3)  => (8.90, 8.24, 0.96),
)

# Cation electronegativity of hydrogen. The quadratic χ⁺ = a + b + c would give 6.85 and let
# hydrogen donate far too much charge; the literature value 20.02 keeps H charges realistic.
const GASTEIGER_IONXH = 20.02

function _gasteiger_mode(sym::Symbol, hyb::Symbol, n_O::Int)
    sym in (:H, :F, :Cl, :Br, :I) && return :any
    if sym == :S && hyb == :none
        return n_O >= 2 ? :so2 : n_O == 1 ? :so : :sp3
    end
    return hyb
end

function _gasteiger_abc(sym::Symbol, mode::Symbol, i::Int)
    haskey(GASTEIGER_PARAMS, (sym, mode)) ||
        throw(ArgumentError("no Gasteiger parameters for atom $i ($sym, $mode); charges are parameterized for common organic elements only"))
    return GASTEIGER_PARAMS[(sym, mode)]
end

"""
    q = partial_charges(mol::SimpleMolGraph; niter=12)

Gasteiger–Marsili (PEOE) partial charge of every atom of `mol`, as a vector over its
vertices. Charges are seeded from the formal charges and iteratively equalized across bonds;
the sum over all atoms (including implicit hydrogens; see below) equals the total formal
charge.

Hydrogens that are graph vertices are charged like any other atom. Hydrogens that are
implicit share one additional charge state per heavy atom, which is folded into that heavy
atom's reported value, so `q[i]` is the charge of atom `i` *plus its implicit hydrogens*.
Use [`atom_group_charges`](@ref) to also fold in explicit hydrogens, giving one charge per
heavy atom.

Throws an `ArgumentError` for elements outside the parameterization (common organic elements
plus halogens, S, and P).
"""
function partial_charges(mol::SimpleMolGraph; niter::Integer=12)
    n = nv(mol.graph)
    syms = atom_symbol(mol)
    hyb = hybridization(mol)
    nimplicit = implicit_hydrogens(mol)

    # oxygen neighbors, for the unperceived-sulfur fallback modes
    n_O = zeros(Int, n)
    for i in 1:n
        syms[i] == :S || continue
        n_O[i] = count(nb -> syms[nb] == :O, neighbors(mol, i))
    end

    abc = [_gasteiger_abc(syms[i], _gasteiger_mode(syms[i], hyb[i], n_O[i]), i) for i in 1:n]
    ionX = [syms[i] == :H ? GASTEIGER_IONXH : sum(abc[i]) for i in 1:n]
    abcH = GASTEIGER_PARAMS[(:H, :any)]

    # Seed with formal charges, shared equally among resonance-equivalent terminal atoms:
    # same-element single-neighbor atoms on a common center (the oxygens of a carboxylate,
    # nitro, sulfonate, or phosphate group) split their total formal charge. This mirrors the
    # conjugated-charge splitting of the RDKit reference implementation for such groups;
    # charges without equivalent partners (e.g. an ammonium nitrogen) stay localized.
    q = Float64.(atom_charge(mol))
    for c in 1:n
        nbrs = neighbors(mol, c)
        for elem in unique(syms[j] for j in nbrs)
            elem == :H && continue
            term = [j for j in nbrs if syms[j] == elem && length(neighbors(mol, j)) == 1]
            length(term) >= 2 || continue
            tot = sum(q[j] for j in term)
            iszero(tot) && continue
            for j in term
                q[j] = tot / length(term)
            end
        end
    end
    hq = zeros(n)                    # charge of one implicit hydrogen of atom i
    χ = zeros(n)
    χh = zeros(n)
    dq = zeros(n)
    dhq = zeros(n)
    damp = 0.5
    for _ in 1:niter
        for i in 1:n
            a, b, c = abc[i]
            χ[i] = a + q[i] * (b + c * q[i])
            if nimplicit[i] > 0
                a, b, c = abcH
                χh[i] = a + hq[i] * (b + c * hq[i])
            end
        end
        fill!(dq, 0.0)
        fill!(dhq, 0.0)
        for e in edges(mol.graph)
            i, j = src(e), dst(e)
            dx = χ[j] - χ[i]
            # charge flows toward the more electronegative atom, scaled by the donor's χ⁺
            t = dx / (dx > 0 ? ionX[i] : ionX[j])
            dq[i] += t
            dq[j] -= t
        end
        for i in 1:n
            nimplicit[i] > 0 || continue
            dx = χh[i] - χ[i]
            t = dx / (dx > 0 ? ionX[i] : GASTEIGER_IONXH)
            dq[i] += nimplicit[i] * t
            dhq[i] -= t
        end
        for i in 1:n
            q[i] += damp * dq[i]
            hq[i] += damp * dhq[i]
        end
        damp *= 0.5
    end

    # report each implicit hydrogen's charge on its heavy atom, so the vector sums to the
    # total formal charge over the graph's vertices
    return [q[i] + nimplicit[i] * hq[i] for i in 1:n]
end

"""
    qg = atom_group_charges(mol::SimpleMolGraph; niter=12)

Per-atom charges from [`partial_charges`](@ref) with every explicit hydrogen's charge folded
onto the heavy atom it is bonded to. Hydrogen vertices report zero, so summing `qg` over
heavy atoms gives the total formal charge. This is the natural amplitude for a per-heavy-atom
Coulomb field.
"""
function atom_group_charges(mol::SimpleMolGraph; niter::Integer=12)
    q = partial_charges(mol; niter)
    syms = atom_symbol(mol)
    qg = copy(q)
    for i in eachindex(qg)
        syms[i] == :H || continue
        heavies = [nb for nb in neighbors(mol, i) if syms[nb] != :H]
        length(heavies) == 1 || continue   # e.g. H₂ or a lone hydride: leave the charge on H
        qg[only(heavies)] += q[i]
        qg[i] = 0.0
    end
    return qg
end
