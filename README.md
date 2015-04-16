# GibbsSeaWater

This is a Julia wrapper for the [Gibbs Seawater Oceanographic Toolbox
(GSW)](http://www.teos-10.org/software.htm#1), which implements the TEOS-10
subroutines for the thermodynamic seawater equation of state.

The current included GSW version is 3.03.

## Installation

Installation should work with

```julia
Pkg.clone(git@github.com:njwilson23/GibbsSeaWater.jl.git)
```

Right now, I suspect that the build script only works on Linux. If someone with
access to OS X and Windows machines could test and make the necessary
modifications to `build.jl`, I'd be happy to accept pull requests.

