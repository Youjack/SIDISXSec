"""
Part of SIDISStructFunc.jl

PV17 model of unpolarized structure functions
"""
module SF_PV17

using TMDlib: mkSF, FxzqTQ
using ..SIDISXSec

export get_sf_pv17

function FUUT(sf::Ptr{Cvoid}, var::SidisVar, μ²::Real)::Float64
    xB, _, _, _, _, _, zh, _, _, PhT² = exposestruct(var)
    qT = √get_qT²(zh, PhT²)
    return FxzqTQ(sf, xB, zh, qT, √μ²)
end

function get_sf_pv17()::SidisStructFunc
    sf = mkSF("PV17_grid_FUUT_Pip/50")
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (data, var, μ², rtol=0.0) -> FUUT(sf, var, μ²),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
    )
end

end # module
