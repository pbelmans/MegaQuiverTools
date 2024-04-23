# MegaQuiverTools

[![tests](https://github.com/Catullo99/MegaQuiverTools.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/Catullo99/MegaQuiverTools.jl/actions/workflows/Runtests.yml)
[![docs](https://github.com/Catullo99/MegaQuiverTools.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/Catullo99/MegaQuiverTools.jl/actions/workflows/Documenter.yml)

This is the stripped down version of QuiverTools that will be submitted to MEGA.

## Installation

Once a functioning version of Julia is installed on the machine, the following command can be used to install the package:

```julia-repl
]add https://github.com/Catullo99/MegaQuiverTools.jl.git
```

One can alternatively clone the repository and install the package from the local copy:

```julia-repl
]dev path/to/MegaQuiverTools.jl
```

This _should_ add the package to the Julia environment and let Julia understand that it is being developed. **Avoid this unless you plan to modify this package and you know what you are doing.**

## Documentation

The documentation is hosted on GitHub Pages. It can be found [here](https://catullo99.github.io/MegaQuiverTools.jl/).

To build a local .pdf version of the documentation, first modify the make.jl file by replacing the parameter ``format = Documenter.HTML()`` with ``format = Documenter.LaTeX()``, then run

```bash
julia --project make.jl
```

from inside the ```docs/``` directory. Note that this requires a working installation of LaTeX.
