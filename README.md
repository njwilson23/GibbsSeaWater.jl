# GibbsSeaWater

<img src="https://travis-ci.org/njwilson23/GibbsSeaWater.jl.svg" alt="Travis CI status">

This is a Julia wrapper for the [Gibbs Seawater Oceanographic Toolbox
(GSW)](http://www.teos-10.org/software.htm#1), which implements the TEOS-10
subroutines for the thermodynamic seawater equation of state.

The current included GSW version is 3.03.

## Installation

Installation should work with

```julia
Pkg.clone(git@github.com:njwilson23/GibbsSeaWater.jl.git)
```

provided that `make` and a C compiler are available.

While I've tried to write this in a way to work on all systems, I've only tested
with Linux. If someone with access to an OS X or Windows machine finds problems,
I'd be happy to accept pull requests.

## Usage

All GSW 3.03 functions are available. For example, to convert between potential
and conservative temperature,

```julia
import GibbsSeaWater
# Sp = 34, pres = 150 dbar, lon = -65, lat = 40
Sa = GibbsSeaWater.sa_from_sp(34.0, 150.0, -65.0, 40.0)  
```

For more information, see the GSW documentation at
[http://www.teos-10.org/pubs/gsw/html/gsw_contents.html](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html).

