"""
Part of SIDISStructFunc.jl

Inverse-Error Weighting [echevarria2018matching]
"""
module InEW

using ..SIDISXSec

export get_sf_InEW

const _rtol = 5e-3

"Inverse-Error Weighting, assuming uniform power corrections to all structure functions."
function FUU_InEW(M, FUU_tmd::Function, FUU_coll::Function, xB, Q², zh, qT², μ², rtol=_rtol, nΔ=0)::Float64
    qT, Q = √qT², √Q²
    W() = FUU_tmd( xB, Q², zh, qT², μ², rtol)
    Z() = FUU_coll(xB, Q², zh, qT², μ², rtol)
    Δ_W = (qT/Q)^2 + (M/Q)^2
    Δ_Z = (M/qT)^2 * ( 1 + 1/4*log(Q²/qT²+1)^2 )
    Δ_Γ = Δ_W * Δ_Z / √( Δ_W^2 + Δ_Z^2 )
    Γ = ( Δ_Z^2 * W() + (qT/Q > 0.01 ? Δ_W^2 * Z() : 0) )/
        ( Δ_W^2 + Δ_Z^2 )
    return ( 1 + nΔ * Δ_Γ ) * Γ
end

"""
    get_sf_InEW(M, sf_tmd::SidisStructFunc, sf_coll::SidisStructFunc, nΔ=0)::SidisStructFunc

Get structure functions with inverse-error weighting.
The error `nΔ * ΔΓ` is added to the functions.
"""
function get_sf_InEW(M, sf_tmd::SidisStructFunc, sf_coll::SidisStructFunc, nΔ=0)::SidisStructFunc
    return SidisStructFunc(
        (xB, Q², zh, qT², μ², rtol=_rtol) -> FUU_InEW(M, sf_tmd.FUUL,      sf_coll.FUUL,      xB, Q², zh, qT², μ², rtol, nΔ),
        (xB, Q², zh, qT², μ², rtol=_rtol) -> FUU_InEW(M, sf_tmd.FUUT,      sf_coll.FUUT,      xB, Q², zh, qT², μ², rtol, nΔ),
        (xB, Q², zh, qT², μ², rtol=_rtol) -> FUU_InEW(M, sf_tmd.FUUcosϕh,  sf_coll.FUUcosϕh,  xB, Q², zh, qT², μ², rtol, nΔ),
        (xB, Q², zh, qT², μ², rtol=_rtol) -> FUU_InEW(M, sf_tmd.FUUcos2ϕh, sf_coll.FUUcos2ϕh, xB, Q², zh, qT², μ², rtol, nΔ),
    )
end

end # module
