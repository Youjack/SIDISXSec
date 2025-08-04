"""
Part of SIDISStructFunc.jl

SV19 model of unpolarized structure functions
"""
module SF_SV19

using SpecialFunctions: besselj
using QuadGK
using ..SIDISXSec
using ..SIDISXSec.QCDData
using TMDTools.SV19

export _get_sf_sv19, get_sf_sv19

const CF = 4/3

const _rtol = 5e-3
# const _qTmin = 1e-6
const _qTmin = 1e-4

#= function FUUT(f̃::Function, D̃::Function, xB, Q², zh, qT², μ², od0::OgataData, rtol=_rtol)::Float64
    qT = max(√qT², _qTmin)
    return xB * 1/(2π)^2 * ft_symmetric2d(od0,
        bT -> let #= bT = √(bT^2 + SV19.b₀²/Q²) =#
            sum(i -> quark_charge[i]^2 *
                f̃(quark_code[i], xB, bT^2, μ², Q², rtol) * D̃(quark_code[i], zh, bT^2, μ², Q², rtol),
                1:num_quark)
        end,
        rtol=rtol)(qT)
end =#
function FUUT(αs::Function, f̃::Function, D̃::Function, xB, Q², zh, qT², μ², rtol=_rtol)::Float64
    qT = √qT²
    return xB *
        ( 1 + CF * αs(μ²)/4π *( - 16 + π^2/3 + 6log(Q²/μ²) - 2log(Q²/μ²)^2 ) ) * # LO + NLO
        1/2π * quadgk(bT -> let #= bT = √(bT^2 + SV19.b₀²/Q²) =#
            bT * besselj(0, qT*bT) * sum(i -> quark_charge[i]^2 *
                f̃(quark_code[i], xB, bT^2, μ², Q², rtol) * D̃(quark_code[i], zh, bT^2, μ², Q², rtol),
                1:num_quark)
        end, 0, Inf, rtol=rtol)[1]
end

"""
    get_sf_sv19(data::SidisData; order=2)::SidisStructFunc

`data` should be the `LHASplitExtend` version.
"""
function get_sf_sv19(data::SidisData; order=2)::SidisStructFunc
    sv19data = SV19.setdata(data.f, data.D, data.αs)
    f̃(q, x, bT², μ², ζ, rtol) = SV19.get_tmdpdf(sv19data, q, μ², ζ, order=order)(x, bT²)
    D̃(q, z, bT², μ², ζ, rtol) = SV19.get_tmdff( sv19data, q, μ², ζ, order=order)(z, bT²)
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=_rtol) -> FUUT(data.αs, f̃, D̃, xB, Q², zh, qT², μ², rtol),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
    )
end

end # module
