# MolecularGaussians.jl

Alignmnent and comparison of small molecules read from .sdf files represented as Gaussian Mixture Models. 

## Build GMMs from molecules

```julia
julia> using MolecularGaussians

julia> using MolecularGraph

julia> # note that by default, all standard deviations are set give Gaussians the same volume as the atom the represent, and all weights are set to 1.0

julia>  mol1=sdftomol(joinpath(dirname(pathof(MolecularGaussians)), "..", "data", "E1050_3d.sdf"));       

julia>  mol2=sdftomol(joinpath(dirname(pathof(MolecularGaussians)), "..", "data", "E1103_3d.sdf"));       

julia> molgmm1 = MolGMM(mol1)
MolGMM from molecule with formula C18H24O8S2 with 52 Gaussians.


julia> molgmm2 = MolGMM(mol2)
MolGMM from molecule with formula C18H24O5S with 48 Gaussians.


julia> fmolgmm1 = FeatureMolGMM(mol1)
FeatureMolGMM from molecule with formula C18H24O8S2 with 13 Gaussians in 4 GMMs with labels:
[:acceptor, :rings, :aromaticrings, :anion]


julia> fmolgmm2 = FeatureMolGMM(mol2)
FeatureMolGMM from molecule with formula C18H24O5S with 10 Gaussians in 5 GMMs with labels:
[:acceptor, :rings, :donor, :anion, :aromaticrings]
```

## Compute L2 distance between two GMMs
```julia
julia> moldist = gmm_distance(molgmm1, molgmm2)
28.463341354279294

julia> fmoldist = gmm_distance(fmolgmm1, fmolgmm2)
10.332058976004422
```

## Find transformation to align GMMs (minimize distance)
```julia
julia> overlap, tform, nevals = tiv_gogma_align_gmms(molgmm1, molgmm2, 2, 2; maxstagnant=100);

julia> gmm_distance(molgmm1, molgmm2; tform2=tform)
15.74638368864396

julia> overlap, tform, nevls = tiv_gogma_align_gmms(fmolgmm1, fmolgmm2);

julia> gmm_distance(fmolgmm1, fmolgmm2, tform2=tform)
6.085914991089279
```