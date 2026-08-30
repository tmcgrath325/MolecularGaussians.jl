```@meta
CurrentModule = MolecularGaussians
```

# API reference

```@docs
MolecularGaussians
```

## Building pharmacophore models

```@docs
PharmacophoreGMM
feature_maps
feature_labels
atoms_to_feature
transform
```

## Bond rotation

```@docs
bondrotate
```

## Flexible alignment

`flexible_align` globally aligns a molecule onto another while rotating about each of its
rotatable bonds, searching the bond angles continuously alongside the rigid transform rather
than enumerating discrete conformers. `aligned` returns the posed, flexed molecule and
`joint_angles` the optimized bond angles.

```@docs
flexible_align
```

## Feature definitions

```@docs
parse_feature_definitions
FamilyDef
FeatureDef
AtomType
smarts_or
smarts_andnot
```

## Comparison and alignment

The comparison and alignment operations are re-exported from
[GaussianMixtureAlignment.jl](https://github.com/tmcgrath325/GaussianMixtureAlignment.jl)
and act on `PharmacophoreGMM`s directly; see that package for their full
documentation, or query them at the REPL with `?`.

- `overlap`, `distance`, `tanimoto` — compare two mixtures (see [Concepts](@ref)).
- `local_align` — refine a pose to a nearby overlap maximum.
- `gogma_align`, `tiv_gogma_align` — global branch-and-bound alignment.
- `flex_gogma_align`, `aligned`, `joint_angles` — the flexible branch-and-bound
  search underlying [`flexible_align`](@ref) and its result accessors.
- `rocs_align` — local alignment from a starting pose.

Because `distance` is also exported by other loaded packages, qualify it as
`MolecularGaussians.distance` (or select any of the exporting modules) to avoid
an ambiguity error.

## Visualization

`moldisplay`, re-exported from
[MolecularGraph.jl](https://github.com/mojaie/MolecularGraph.jl), draws a
molecule. `pharmacophoredisplay` overlays one or more `PharmacophoreGMM`s on
their molecules and requires a Makie backend.

```@docs
pharmacophoredisplay
```
