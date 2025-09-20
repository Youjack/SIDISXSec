"""
Part of SIDISStructFunc.jl

Collinear factorization
"""
module Coll

using QuadGK
using LinearAlgebra: ⋅
using ..SIDISXSec
using ..SIDISXSec.Constants, ..SIDISXSec.QCDData
using Printf

export _get_sf_coll, _get_sf_coll_nlp
export get_sf_coll

const CF = 4/3
const TF = 1/2

const _rtol = 5e-3
const _sqrteps = √eps(Float64)

#= (real-time calculation) ========================================================================#

# C = ( C_q_qg, C_q_gq, C_g_qq̄ )
CUUT(x̂, ẑ) = (
    2CF *( (1-x̂) * (1-ẑ) + ( 1 + x̂^2 *    ẑ ^2 )/( (1-x̂) * (1-ẑ) ) ),
    2CF *( (1-x̂) *    ẑ  + ( 1 + x̂^2 * (1-ẑ)^2 )/( (1-x̂) *    ẑ  ) ),
    2TF * 1/(ẑ*(1-ẑ)) * ( x̂^2 + (1-x̂)^2 ) * ( ẑ^2 + (1-ẑ)^2 )       ,
)
CUUcosϕh(x̂, ẑ) = (
    -4CF * ( x̂ *    ẑ  + (1-x̂) * (1-ẑ) ) *           √( x̂ * ẑ /( (1-x̂) * (1-ẑ) ) ),
     4CF * ( x̂ * (1-ẑ) + (1-x̂) *    ẑ  ) * (1-ẑ)/ẑ * √( x̂ * ẑ /( (1-x̂) * (1-ẑ) ) ),
    -4TF * (2x̂-1) * (2ẑ-1)               * (1-x̂)/ẑ * √( x̂ * ẑ /( (1-x̂) * (1-ẑ) ) ), # can be positive/negative
)
CUUcos2ϕh(x̂, ẑ) = (
    4CF * x̂ *    ẑ ,
    4CF * x̂ * (1-ẑ),
    8TF * x̂ * (1-x̂),
)
CUUL(x̂, ẑ) = 2 .* CUUcos2ϕh(x̂, ẑ)
CLL(x̂, ẑ) = (
    2CF *( 2(x̂+ẑ) + ( x̂^2 + ẑ^2 )/( (1-x̂) * (1-ẑ) ) ),
    2CF *( 2x̂ + 2(1-ẑ) + ( x̂^2 + (1-ẑ)^2 )/( (1-x̂) * ẑ ) ),
    2TF * (2x̂-1) * ( ẑ^2 + (1-ẑ)^2 )/( ẑ * (1-ẑ) ),
)

function _C_f_D(αs::Function, f::Function, D::Function, CUU::Function,
        xB, zh, qT²oQ², μ², rtol)::Float64
    zmin = zh *( 1 + xB/(1-xB) * qT²oQ² )
    if zmin ≥ 1.0 return 0.0 end
    return 4π * αs(μ²) * quadgk(
        z -> let ẑ = zh/z, x = xB *( 1 + ẑ/(1-ẑ) * qT²oQ² ), x̂ = xB/x
            if (x < 1 + _sqrteps) x = min(x, 1.0) end
            1/z * (x̂*ẑ)/(1-ẑ) *
            sum(i -> quark_charge[i]^2 *( CUU(x̂,ẑ) ⋅ ( # ↓ force positive ↓
                max(1e-10, f(quark_code[i], x, μ²)) * max(1e-10, D(quark_code[i], z, μ²)),
                max(1e-10, f(quark_code[i], x, μ²)) * max(1e-10, D(21           , z, μ²)),
                max(1e-10, f(21           , x, μ²)) * max(1e-10, D(quark_code[i], z, μ²)),
                )), 1:num_quark)
        end,
        zmin, 1.0, rtol=rtol)[1]
end
function _F_coll(αs::Function, f::Function, D::Function, C::Function,
        xB, Q², zh, qT², μ², rtol=_rtol)::Float64
    1/(2π)^3 * xB/(2zh^2*Q²) * _C_f_D(αs, f, D, C, xB, zh, qT²/Q², μ², rtol)
end

function _get_sf_coll(data::SidisData)::SidisStructFunc
    αs, f, g, D = data.αs, data.f, data.g, data.D
    return SidisStructFunc(
        (xB, Q², zh, qT², μ², rtol=_rtol) -> _F_coll(αs, f, D, CUUL,      xB, Q², zh, qT², μ², rtol),
        (xB, Q², zh, qT², μ², rtol=_rtol) -> _F_coll(αs, f, D, CUUT,      xB, Q², zh, qT², μ², rtol),
        (xB, Q², zh, qT², μ², rtol=_rtol) -> _F_coll(αs, f, D, CUUcosϕh,  xB, Q², zh, qT², μ², rtol),
        (xB, Q², zh, qT², μ², rtol=_rtol) -> _F_coll(αs, f, D, CUUcos2ϕh, xB, Q², zh, qT², μ², rtol),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=_rtol) -> _F_coll(αs, g, D, CLL,       xB, Q², zh, qT², μ², rtol),
    )
end

#= (collgrid) =====================================================================================#

include("../TMDGrid/generate_collgrid.jl")

function F_coll(F_grid::Function, xB, Q², zh, qT², μ²)::Float64
    qT²oQ², R, A = grid_encode(xB, zh, qT²/Q²)
    round3(x) = round(x, sigdigits=3)
    # use ≈ ?
    # μ²
    if μ² < μ²_low - _sqrteps || μ² > μ²_up + _sqrteps
        throw(ErrorException("μ² = $(round3(μ²)) out of bounds [$(round3(μ²_low)),$(round3(μ²_up))]."))
    else
        μ² = clamp(μ², μ²_low, μ²_up)
    end
    # qT²/Q²
    if qT²oQ² < qT²oQ²_low - _sqrteps || qT²oQ² > qT²oQ²_up′ + _sqrteps
        throw(ErrorException("qT²/Q² = $(round3(qT²/Q²)) out of bounds [$(round3(qT²oQ²_low)),$(round3(qT²oQ²_up′))]."))
    else
        qT²oQ² = clamp(qT²oQ², qT²oQ²_low, qT²oQ²_up′)
    end
    # xB, zh
    if xB < xB_low - _sqrteps || xB > 1 + _sqrteps
        throw(ErrorException("xB = $(round3(xB)) out of bounds [$(round3(xB_low)),1.0]."))
    else
        xB = clamp(xB, xB_low, 1)
    end
    if zh < zh_low - _sqrteps || zh > 1 + _sqrteps
        throw(ErrorException("zh = $(round3(zh)) out of bounds [$(round3(zh_low)),1.0]."))
    else
        zh = clamp(zh, zh_low, 1)
    end
    qT²oQ²_max = get_qT²oQ²max(xB,zh)
    if qT²oQ² ≥ qT²oQ²_max
        return 0
    end
    # return
    R = clamp(R, 0,1)
    A = qT²oQ²_max == qT²oQ²_up ? 0.0 : clamp(A, 0,1)
    return F_grid(qT²/Q², R, A, μ²)
end

function get_sf_coll(name::String)::SidisStructFunc
    FUUL_grid     = interpolate_tmdgrid(read_tmdgrid("Coll_FUUL_$name"))
    FUUT_grid     = interpolate_tmdgrid(read_tmdgrid("Coll_FUUT_$name"))
    FUUcosϕh_grid = interpolate_tmdgrid(read_tmdgrid("Coll_FUUcosϕh_$name"))
    FLL_grid      = interpolate_tmdgrid(read_tmdgrid("Coll_FLL_$name"))
    return SidisStructFunc(
        (xB, Q², zh, qT², μ², rtol=0.0) -> F_coll(FUUL_grid,     xB, Q², zh, qT², #= μ² =#Q²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> F_coll(FUUT_grid,     xB, Q², zh, qT², #= μ² =#Q²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> F_coll(FUUcosϕh_grid, xB, Q², zh, qT², #= μ² =#Q²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> F_coll(FUUL_grid,     xB, Q², zh, qT², #= μ² =#Q²)/2,
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=0.0) -> F_coll(FLL_grid,      xB, Q², zh, qT², Q²),
    )
end

end # module
