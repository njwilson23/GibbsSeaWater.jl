# GibbsSeaWater

<a href="https://travis-ci.org/njwilson23/GibbsSeaWater.jl">
    <img src="https://travis-ci.org/njwilson23/GibbsSeaWater.jl.svg" alt="Travis CI status">
</a>

This is a Julia wrapper for the [Gibbs Seawater Oceanographic Toolbox
(GSW)](http://www.teos-10.org/software.htm#1), which implements the TEOS-10
subroutines for the thermodynamic seawater equation of state.

The current included GSW version is 3.03.

## Installation

Installation should work with

```julia
Pkg.add("GibbsSeaWater")
```

provided that `make` and a C compiler are available.

While I've tried to write this in a way to work on all systems, I've only tested
with Linux. If someone with access to an OS X or Windows machine finds problems,
I'd be happy to accept pull requests.

## Usage

All GSW 3.03 functions are available. For example, to convert between practical
salinity (PSS-78) and absolute salinity (TEOS-10),

```julia
import GibbsSeaWater
# Sp = 34, pres = 150 dbar, lon = -65, lat = 40
Sa = GibbsSeaWater.sa_from_sp(34.0, 150.0, -65.0, 40.0)  
```

For more information, see the GSW documentation at
[http://www.teos-10.org/pubs/gsw/html/gsw_contents.html](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html).

