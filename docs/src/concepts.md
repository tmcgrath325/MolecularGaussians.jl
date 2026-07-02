```@meta
CurrentModule = MolecularGaussians
```

# Concepts

This page explains the representation MolecularGaussians.jl uses and the pieces
that build it up.

## Molecules as Gaussian mixtures

A molecule is a set of atoms with 3-D coordinates. MolecularGaussians.jl places
an isotropic Gaussian on each chosen point of interest — a heavy atom, or the
centroid of a group of atoms — giving each a width `σ` and an amplitude `ϕ`. The
collection of Gaussians is an `IsotropicGMM` (from
[GaussianMixtureAlignment.jl](https://github.com/tmcgrath325/GaussianMixtureAlignment.jl)),
and comparing two molecules reduces to comparing two smooth density functions
rather than matching discrete atoms.

Two widths are available. A *volume* Gaussian approximates an atom's van der
Waals sphere, so the mixture models molecular shape; this is the default when a
[`PharmacophoreGMM`](@ref) is built with no feature map. A *feature* Gaussian
instead summarizes a pharmacophore group. [`atoms_to_feature`](@ref) builds a
single Gaussian from one index set: a lone atom takes its own amplitude and
width, while a multi-atom set is centered at the atoms' centroid with a width
derived from their combined van der Waals volume.

## Pharmacophore features and families

A *pharmacophore feature* is a chemical group responsible for a type of
interaction — a hydrogen-bond donor or acceptor, an aromatic ring, an ionizable
group, and so on. Each such type is a *family*, named by a `Symbol` such as
`:Donor` or `:Acceptor`.

Families are specified by SMARTS substructure patterns. The building blocks are:

- [`AtomType`](@ref): a named atom class, given by the SMARTS expression inside a
  bracketed atom. [`smarts_or`](@ref) and [`smarts_andnot`](@ref) combine atom
  types into more specific classes.
- [`FeatureDef`](@ref): a SMARTS pattern together with the family it defines and
  a per-atom weight vector.
- [`FamilyDef`](@ref): a complete parsed set of atom types, feature definitions,
  and the grouping of feature names into families.

[`parse_feature_definitions`](@ref) reads a `.fdef` file — by default the
definitions bundled with the package — and returns a `FamilyDef`.

## From molecule to `PharmacophoreGMM`

Building a pharmacophore model proceeds in two stages:

1. [`feature_maps`](@ref) matches the SMARTS patterns of a `FamilyDef` against a
   molecule and returns, for each requested family, the sets of atom indices
   that realize it. Distinct matches over the same atoms are collapsed to one
   set.
2. [`PharmacophoreGMM`](@ref) turns each index set into a single Gaussian via
   [`atoms_to_feature`](@ref), tagging each with its family. The result stores
   the Gaussians and their parallel family labels alongside the underlying
   molecular graph.

## Comparison and alignment

Because every Gaussian in a `PharmacophoreGMM` carries a family label, overlap
is restricted to Gaussian pairs whose labels interact — by default only pairs
sharing a family, so comparisons are made family-by-family and summed:

- `overlap` integrates the product of two mixtures — larger when their features
  coincide in space and type.
- `distance` is the `L₂` distance between the two density functions.
- `tanimoto` normalizes overlap into a similarity score in `[0, 1]`.

These quantities depend on the relative pose of the two molecules. Alignment
searches for the rigid transformation (rotation and translation) that maximizes
overlap:

- `gogma_align` and `tiv_gogma_align` perform a global branch-and-bound search.
- `rocs_align` performs a local optimization from a starting pose.
- `local_align` refines a pose to a nearby overlap maximum.

Applying an alignment result's transformation to a `PharmacophoreGMM` — with the
result's `tform` or with [`transform`](@ref) — moves both its Gaussians and its
molecular graph, so the aligned model can be compared or plotted directly.
```
