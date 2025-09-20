"""
Part of SIDISStructFunc.jl

SV19 model of unpolarized structure functions
"""
module SF_SV19

using SpecialFunctions: besselj
using QuadGK
using ..SIDISXSec
using ..SIDISXSec.QCDData
using TMDTools.SV19, TMDTools.ZetaEvolve
using ..SIDISXSec.TMD: HUUT

export get_sf_sv19

const _rtol = 5e-3

function FUUT(f̃::Function, D̃::Function, xB, Q², zh, qT², μ², rtol=_rtol)::Float64
    qT = √qT²
    return xB *
        1/2π * quadgk(bT -> let
            bT * besselj(0, qT*bT) * sum(i -> quark_charge[i]^2 *
                f̃(quark_code[i], xB, bT^2, μ², Q², rtol) * D̃(quark_code[i], zh, bT^2, μ², Q², rtol),
                1:num_quark)
        end, 0, Inf, rtol=rtol)[1]
end

"""
    get_sf_sv19(data::SidisData; order=2, inclH=false)::SidisStructFunc

- `data` should be the `LHASplitExtend` version.
- set `inclH=true` to include hard parts.
"""
function get_sf_sv19(data::SidisData; order=2, inclH=false)::SidisStructFunc
    sv19data = SV19.setdata(data.f, data.D, data.αs)
    f̃(q, x, bT², μ², ζ, rtol) = SV19.get_tmdpdf(sv19data, q, μ², ζ, order=order)(x, bT²)
    D̃(q, z, bT², μ², ζ, rtol) = SV19.get_tmdff( sv19data, q, μ², ζ, order=order)(z, bT²)
    _HUUT(Q², μ²) = inclH ? HUUT(data.αs, Q², μ²) : 1.
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=_rtol) -> _HUUT(Q², μ²) * FUUT(f̃, D̃, xB, Q², zh, qT², μ², rtol),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
    )
end

end # module
