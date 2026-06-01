# FTAN.jl

Array-first Julia tools for Multiple Filter Technique (MFT/FTAN) analysis of
surface-wave waveforms.

This package is split out from the working `MFT_v2.jl` Pluto notebook. The
package API intentionally uses arrays plus explicit `dt`, `distance`, and
`periods`; it does not include the old `SeismicTrace` wrapper API.

## Local use

```julia
import Pkg
Pkg.develop(path="/path/to/FTAN.jl")
using FTAN
```

## Documentation notebooks

Runnable Pluto examples live in `pluto_notebooks/`. They are designed to export
as static documentation through GitHub Pages using PlutoSliderServer.

The GitHub Pages workflow uses the current Julia 1.12 stable line and installs
the current registered Pluto/PlutoSliderServer releases during export.
