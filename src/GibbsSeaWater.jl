module GibbsSeaWater

# package code goes here
const _MODULEPATH = Pkg.dir("GibbsSeaWater")
const _LIBGSWPATH = joinpath(MODULEPATH, "deps/build/libgswteos-10")

function sigma0(sa, ct)
    return ccall((:gsw_sigma0, _LIBGSWPATH), Float64, (Float64, Float64),
                 sa, ct)
end

end # module
