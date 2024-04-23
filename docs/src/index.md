# MegaQuiverTools.jl

## Contents

```@contents
```

## Introduction

MegaQuiverTools is a beta version of QuiverTools,
a software suite in development for treatment of quiver representations,
their roots, their moduli spaces and computations of several invariants of these.

Examples of the functionalities of MegaQuiverTools will be:

- Treatment of Schur roots, canonical decomposition of dimension vectors;
- Computation of the Harder-Narasimhan stratification of the parameter space of quiver moduli;
    - Computation of weights of various vector bundles relative to the 1-parameter subgroups associated to each stratum;
    - Verification of Teleman quantization criteria;
- Construction of the Chow ring of a given quiver moduli, computation of Euler characteristics;
- Computation of the Hodge diamond of quiver moduli;
- etc.

QuiverTools will be available as a Julia package and as a Sage library.

Features present in this beta version are

 - Treatment of Harder-Narasimhan stratifications, weight computations;
 - Chow rings and Euler characteristic computations.

## Installation

At the moment the package is not registered, so you can install it by running the following command in the Julia REPL:

```julia-repl
pkg>add https://github.com/Catullo99/MegaQuiverTools.jl.git
```


## Aknowledgements

QuiverTools is being developed by [P. Belmans](https://pbelmans.ncag.info/), H. Franzen and [G. Petrella](https://www.giannipetrella.eu).

The Julia version is developed and maintained by [G. Petrella](https://www.giannipetrella.eu).

P. B. was partially supported by the Luxembourg National Research Fund (FNR–17113194).

H. F. was partially supported by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) SFB-TRR 358/1 2023 “Integral Structures in Geometry and Representation Theory” (491392403).

G. P. was supported by the Luxembourg National Research Fund (FNR–17953441).