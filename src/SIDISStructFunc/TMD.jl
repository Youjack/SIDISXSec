"""
Part of SIDISStructFunc.jl

Interface for grids of TMD factorized structure functions
"""
module TMD

using QuadGK
using ..SIDISXSec
using ..SIDISXSec.Constants, ..SIDISXSec.QCDData
using Printf

export get_sf_tmd

#= Hard parts =====================================================================================#

const CF = 4/3

# `αs(μ²)`
HUUT(αs, Q², μ²) = 1 + CF * αs(μ²)/4π *( - 16 + π^2/3 + 6log(Q²/μ²) - 2log(Q²/μ²)^2 )

#= (tmdgrid) ======================================================================================#

include("../TMDGrid/generate_tmdgrid.jl")

function get_sf_tmd(name::String)::SidisStructFunc
    FUUT_grid = interpolate_tmdgrid(read_tmdgrid("FUUT_TMD_$name"))
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUT_grid(max(0.0, grid_encode(Q²,qT²)[1]), xB, zh, Q²),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
    )
end
"""
    get_sf_tmd(name::String, αs::Function)::SidisStructFunc

Get the TMD structure functions from TMDGrid `name`.
Pass `αs(μ²)` to get results with hard parts.
"""
function get_sf_tmd(name::String, αs::Function)::SidisStructFunc
    sf = get_sf_tmd(name)
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=0.0) -> HUUT(αs, Q², #= μ² =#Q²) * sf.FUUT(xB, Q², zh, qT², μ²),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
    )
end

end # module
