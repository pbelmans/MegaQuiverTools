using Documenter
using MegaQuiverTools

# memo for myself: when changing dependencies of MegaQuiverTools,
# update the GLOBAL Julia environment to have the documenter work.

```@meta
CurrentModule = MegaQuiverTools
```

makedocs(
    sitename = "MegaQuiverTools",
    format = Documenter.HTML(),
    # format = Documenter.LaTeX(), # builds pdf, does not like the github Documenter action for now. Use only in local build.
    modules = [MegaQuiverTools],
    pages = [   "MegaQuiverTools" => "index.md", 
                "Tutorial" => "tutorial.md",
                "All methods" => "methods.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

# deploydocs(
    # repo = "github.com/Catullo99/MegaQuiverTools.jl.git"
# )