## Physics-based fields: per-atom Gaussian feature families approximating physical
## interactions, as an alternative to SMARTS-derived pharmacophore families. Each heavy atom
## carries one slot per field, so two molecules are compared field-by-field through the
## stacked-GMM overlap and its `interactions` weights:
##
## - `:VdW` — an attractive volume field sized like the default `:Volume` features, rewarding
##   shape complementarity at van der Waals contact.
## - `:Steric` — a narrower Gaussian meant to be weighted *negatively* between models,
##   penalizing interpenetration; together with `:VdW` this is a two-Gaussian caricature of a
##   Lennard-Jones profile.
## - `:Charge` — amplitude equal to the atom's Gasteiger group charge (heavy atom plus its
##   hydrogens; see [`atom_group_charges`](@ref)), so overlap weights carry the sign of
##   q₁·q₂. With a positive interaction coefficient the overlap rewards *matching* charge
##   (electrostatic similarity, as in ligand overlay); use a negative coefficient for
##   complementarity instead.

"""
    physics_feature_maps(mol; families=(:VdW, :Steric, :Charge))

Feature map placing one single-atom feature of every requested physics family on each heavy
atom of `mol`, for use as the `feature_maps` argument of [`PharmacophoreGMM`](@ref).
"""
physics_feature_maps(mol::SimpleMolGraph; families=(:VdW, :Steric, :Charge)) =
    Dict{Symbol,Vector{Vector{Int}}}(f => [[i] for i in heavy_atom_idxs(mol)] for f in families)

"""
    physics_sigma_functions(; steric_scale=0.5, charge_sigma=1.0)

Per-family width functions for the physics fields: `:VdW` uses the volume convention
(`vdw_volume_sigma`), `:Steric` is `steric_scale` times the van der Waals radius, and
`:Charge` has the fixed width `charge_sigma` (in the coordinate units, Å for SDF input).
The defaults are heuristic starting points, not fitted values.
"""
physics_sigma_functions(; steric_scale=0.5, charge_sigma=1.0) = Dict{Symbol,Any}(
    :VdW => vdw_volume_sigma,
    :Steric => (atom, ϕ) -> steric_scale * vdw_radius(atom),
    :Charge => (atom, ϕ) -> charge_sigma,
)

"""
    physics_phi_functions(mol; amplitude=2.7, niter=12)

Per-family amplitude functions for the physics fields of `mol`: constant `amplitude` for
`:VdW` and `:Steric`, and the signed Gasteiger group charge for `:Charge`. The `:Charge`
function closes over charges computed for `mol`, so the returned dictionary must only be
used to build a model of that molecule.
"""
function physics_phi_functions(mol::SimpleMolGraph; amplitude=rocs_amplitude, niter::Integer=12)
    qg = atom_group_charges(mol; niter)
    return Dict{Symbol,Any}(
        :VdW => a -> amplitude,
        :Steric => a -> amplitude,
        :Charge => AtomIndexed((m, i) -> qg[i]),
    )
end

"""
    physics_interactions(; vdw=1.0, steric=1.0, charge=1.0)

Interaction coefficients pairing each physics family with itself, for the `interactions`
keyword of `overlap` and the alignment functions: `vdw` (attractive), `-steric` (the steric
field repels), and `charge` (positive rewards matching charge sign — electrostatic
similarity; pass a negative value for complementarity).
"""
physics_interactions(; vdw=1.0, steric=1.0, charge=1.0) = Dict{Tuple{Symbol,Symbol},Float64}(
    (:VdW, :VdW) => vdw,
    (:Steric, :Steric) => -steric,
    (:Charge, :Charge) => charge,
)

"""
    pgmm = physics_gmm(mol; families=(:VdW, :Steric, :Charge), rigid=false,
                            steric_scale=0.5, charge_sigma=1.0, amplitude=2.7, niter=12)

Build a [`PharmacophoreGMM`](@ref) of `mol` whose features are physics fields — one slot per
requested family on every heavy atom — instead of SMARTS pharmacophore families. Score or
align two such models with the coefficients from [`physics_interactions`](@ref):

    overlap(physics_gmm(x), physics_gmm(y); interactions=physics_interactions())

Signed `:Charge` amplitudes are not supported by the TIV alignment paths (`tiv_gogma_align`
scores with `√(ϕx·ϕy)`); use `gogma_align`, `flexible_align`, or `local_align`.
"""
function physics_gmm(mol::SimpleMolGraph; families=(:VdW, :Steric, :Charge), rigid=false,
                     steric_scale=0.5, charge_sigma=1.0, amplitude=rocs_amplitude, niter::Integer=12)
    σfun = physics_sigma_functions(; steric_scale, charge_sigma)
    ϕfun = :Charge in families ? physics_phi_functions(mol; amplitude, niter) :
        Dict{Symbol,Any}(f => (a -> amplitude) for f in families)
    return PharmacophoreGMM(mol; rigid, feature_maps=physics_feature_maps(mol; families), σfun, ϕfun)
end

# `public` is a soft keyword only on Julia 1.11+; guard so the module still
# loads on the 1.10 LTS, where the declaration is simply absent.
@static if VERSION >= v"1.11"
    eval(Meta.parse("public physics_feature_maps, physics_sigma_functions, physics_phi_functions, atom_group_charges"))
end
