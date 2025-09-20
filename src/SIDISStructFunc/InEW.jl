"""
Part of SIDISStructFunc.jl

Inverse-Error Weighting [echevarria2018matching]
"""
module InEW

using ..SIDISXSec

export get_sf_InEW

const _rtol = 5e-3

"Inverse-Error Weighting, assuming uniform power corrections to all structure functions."
function F_InEW(M, F_tmd::Function, F_coll::Function, xB, Q², zh, qT², μ², rtol=_rtol,
        nΔ=0, Wpow=1, Wnorm=1, Wcut=1)::Float64
    qT, Q = √qT², √Q²
    W() = Wnorm * F_tmd( xB, Q², zh, qT², μ², rtol)
    Z() = F_coll(xB, Q², zh, qT², μ², rtol)
    Δ_W = (Wpow * qT/Q)^2 + (M/Q)^2
    Δ_Z = (M/qT)^2 * ( 1 + 1/4*log(Q²/qT²+1)^2 )
    Δ_Γ = Δ_W * Δ_Z / √( Δ_W^2 + Δ_Z^2 )
    Γ = ( (qT < Wcut*Q ? Δ_Z^2 * W() : 0)
        + (qT/Q > 0.01 ? Δ_W^2 * Z() : 0) )/
        ( Δ_W^2 + Δ_Z^2 )
    return ( 1 + nΔ * Δ_Γ ) * Γ
end

"""
    get_sf_InEW(M, sf_tmd::SidisStructFunc, sf_coll::SidisStructFunc;
        nΔ=0, λ=1, κ=1)::SidisStructFuncs

Get structure functions with inverse-error weighting.
- The error `nΔ * ΔΓ` is added to the functions.
- `Wpow` is a factor to tune the size of power corrections to the TMD term.
- `Wnorm` is a factor to tune the normalization of the TMD term.
- `Wcut` is a factor to tune the cutoff in `qT` of the TMD term.
"""
function get_sf_InEW(M, sf_tmd::SidisStructFunc, sf_coll::SidisStructFunc;
        nΔ=0, Wpow=1, Wnorm=1, Wcut=1)::SidisStructFunc
    InEWargs = (nΔ, Wpow, Wnorm, Wcut)
    return SidisStructFunc(
        (Fargs...) -> F_InEW(M, sf_tmd.FUUL,      sf_coll.FUUL,      Fargs..., InEWargs...),
        (Fargs...) -> F_InEW(M, sf_tmd.FUUT,      sf_coll.FUUT,      Fargs..., InEWargs...),
        (Fargs...) -> F_InEW(M, sf_tmd.FUUcosϕh,  sf_coll.FUUcosϕh,  Fargs..., InEWargs...),
        (Fargs...) -> F_InEW(M, sf_tmd.FUUcos2ϕh, sf_coll.FUUcos2ϕh, Fargs..., InEWargs...),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
        (Fargs...) -> F_InEW(M, sf_tmd.FLL,       sf_coll.FLL,       Fargs..., InEWargs...),
    )
end

end # module
