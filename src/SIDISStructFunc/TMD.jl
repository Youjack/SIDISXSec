"""
Part of SIDISStructFunc.jl

Interface for grids of TMD factorized structure functions
"""
module TMD

using QuadGK
using ..SIDISXSec
using ..SIDISXSec.Constants, ..SIDISXSec.QCDData
using Printf
using SpecialFunctions: zeta

export get_sf_tmd

#= Hard parts =====================================================================================#

const CF = 4/3
const CA = 3

const ζ3 = zeta(3)

# copied from TMDTools.jl
const nfmax = 6
const thresholds = ( 0.00216^2, 0.00467^2, 0.0934^2, 1.27^2, 4.18^2, 172.69^2 )
function get_nf(μ²)
    for i = nfmax:-1:3
        if μ² > thresholds[i] return i end
    end
    throw(DomainError(μ², "get_nf with μ smaller than strange threshold."))
end

# `αs(μ²)`
HUUT(αs, Q², μ², Horder) = let
    H::Float64 = 1
    if Horder > 1
        a = αs(μ²)/4π
        L = log(Q²/μ²)
        H += a * CF *( - 16 + π^2/3 + 6L - 2L^2 )
    end
    if Horder > 2
        H += a^2 *(
            + CF^2 * ( 511/4 + 13*pi^2/3 - 13*pi^4/30 - 60*ζ3 + L * (-93 - 2*pi^2 + 48*ζ3) + L^2 * (-2*pi^2/3 + 50) - 12*L^3 + 2*L^4 )
            + CA * CF * ( -51157/324 - 337*pi^2/54 + 22*pi^4/45 + 626*ζ3/9 + L * (2545/27 + 22*pi^2/9 - 52*ζ3) + L^2 * (2*pi^2/3 - 233/9) + 22*T^3/9 )
            + CF * get_nf(μ²) * ( 4085/162 + 23*pi^2/27 + 4*ζ3/9 - L * (418/27 + 4*pi^2/9) + 38*L^2/9 - 4*L^3/9 )
        )
    end
    return H
end

#= (tmdgrid) ======================================================================================#

include("../TMDGrid/generate_tmdgrid.jl")

"""
    get_sf_tmd(name::String[, αs::Function; Horder=1])::SidisStructFunc

Get the TMD structure functions from TMDGrid `name`.
Pass `αs(μ²)` to get results with the hard factor of order `αs^Horder`.
"""
function get_sf_tmd(name::String, αs::Function; Horder=1)::SidisStructFunc
    FUUT_grid = interpolate_tmdgrid(read_tmdgrid("TMD_FUUT_$name"))
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=0.0) -> HUUT(αs, Q², #= μ² =#Q², Horder) * FUUT_grid(max(0.0, grid_encode(Q²,qT²)[1]), xB, zh, Q²),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
    )
end
function get_sf_tmd(name::String)::SidisStructFunc
    return get_sf_tmd(name, zero, Horder=0)
end

end # module
