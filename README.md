# XCrySDenStructureFormat.jl

This package provides read / write functionality for [XCrySDen XSF](http://www.xcrysden.org/doc/XSF.html) atomic structure files.

It is **strongly** recommended **not** to use this package directly, but rather through [AtomsIO.jl](https://github.com/mfherbst/AtomsIO.jl), which provides a uniform interface (based on [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl)) for reading and writing a large range of atomistic structure files.

## Feature support

Currently supports

- r/w of molecular structures (`ATOMS`)
- r/w of periodic structures (`POLYMER` 1D, `SLAB` 2D, `CRYSTAL` 3D)
- r/w of molecular/periodic trajectories (`.axsf` / `ANIMSTEPS`)
- r/w of forces, using the data key `:force` in parsed `AtomsBase.Atom` instances

Currently does _not_ support

- Data grids (2D and 3D)
- Band grids (`.bxsf`)


## Installation

This package is registered in the General registry, so installation of the latest stable release is as simple as pressing `]` to enter `pkg>` mode in the Julia REPL, and then entering:

```julia
pkg> add XCrySDenStructureFormat
```

or for the development version:

```julia
pkg> dev https://github.com/azadoks/XCrySDenStructureFormat.jl
```
