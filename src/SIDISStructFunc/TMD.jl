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

end # module
