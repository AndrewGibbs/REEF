# REEF (*R*esidue *E*nhanced *E*mbedding *F*ormulae)

# Background
For certain scattering geometries, embedding formulae [1] can be used to efficiently compute the far-field pattern induced by a large range of incident angles.
In principle, only a relatively small subset of far-field patterns, induced by _canonical_ angles are required as inputs to each embedding formula.
Despite the embedding formulae being exact in theory, in [2] it was shown that they are very sensitive to numerical errors in the canonical far-fields.
This Matlab package is a numerically stable implementation of the embedding formulae of [1].
Numerical stability is achieved by adding residue contributions to the embedding formulae, hence the name _Residue Enhanced Embedding Formulae_.

# Basic usage
The user must provide ```M```, and ```p```, which are defined as in [1], the wavenumber ```kwave```, and M canonical far-field patterns induced by incident angles in a vector ```alphas```.
The far-field patterns must be defined as cell array of function handles:

```val = D{m}[obs]```

where ```m``` is the index of the far-field induced by ```alphas(m)```; ```obs``` and ```vals``` are Nx1 vectors of observation angles. Then the code

```
E = Reef(D,alphas,kwave,p);
Eout = E.getFarField(obs_test,inc_test);
```
will efficiently  compute the cross-section for a large number of incident angles ```inc_test``` at the observation angles ```obs_test```.

A full example is provided in EG1.m, which produces the following cross-section with wavenumber 15 on a unit square:
![HNABEMLAB](https://github.com/AndrewGibbs/REEF/blob/main/egsquare.png?raw=true)

# Bibliography

[1] [**Biggs, N.R.T., 2006. A new family of embedding formulae for diffraction by wedges and polygons. Wave Motion, 43(7), pp.517-528.**](https://www.sciencedirect.com/science/article/abs/pii/S0165212506000229)<br>
[2] [**Gibbs, A., Langdon, S. and Moiola, A., 2018. Numerically stable computation of embedding formulae for scattering by polygons. arXiv preprint arXiv:1805.08988.**](https://arxiv.org/abs/1805.08988)
