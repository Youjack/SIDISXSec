"""
Part of SIDISStructFunc.jl

W+Y construction [collins2016relating]
"""
module WY

using ..SIDISXSec
using ..SIDISXSec.QCDData
using TMDTools.MellinConv
using SpecialFunctions: digamma, besselk

export get_sf_WY

const CF = 4/3
const TF = 1/2

const _rtol = 5e-3

#= (real-time calculation) ========================================================================#
# [echevarria2016unpolarized]

#= Collinear to TMD matching coefficients =#
# C = ∑ₙ₌₀ (α/2π)ⁿ Cₙ
MCqq₁(q, μ, ζ)   = CF * 1/(π*q^2) * CoeffFunc(zero, x -> (1+x^2)/(1-x), - log(q^2/ζ) - 3/2)
MCqg₁(q, μ, ζ)   = TF * 1/(π*q^2) * CoeffFunc(x -> (x^2 + (1-x)^2))
z²Mℂqq₁(q, μ, ζ) = CF * 1/(π*q^2) * CoeffFunc(zero, z -> (1+z^2)/(1-z), - log(q^2/ζ) - 3/2)
z²Mℂgq₁(q, μ, ζ) = CF * 1/(π*q^2) * CoeffFunc(z -> (1 + (1-z)^2)/z)

#= Modified collinear to TMD matching coefficients in [collins2016relating] with C₅=1 =#
# const b₀ = 2exp(digamma(1))
# MCqq₁(q, μ, ζ)   = CF/π * CoeffFunc(zero, x -> (1+x^2)/(1-x)   * b₀/(μ*q) * besselk(1,b₀*q/μ), b₀/(μ*q) * besselk(1,b₀*q/μ) *( - log(μ^2/ζ) + log(μ/q) - 3/2 ) + 1/q^2 * besselk(0,b₀*q/μ))
# MCqg₁(q, μ, ζ)   = TF/π * CoeffFunc(      x -> (x^2 + (1-x)^2) * b₀/(μ*q) * besselk(1,b₀*q/μ))
# z²Mℂqq₁(q, μ, ζ) = CF/π * CoeffFunc(zero, z -> (1+z^2)/(1-z)   * b₀/(μ*q) * besselk(1,b₀*q/μ), b₀/(μ*q) * besselk(1,b₀*q/μ) *( - log(μ^2/ζ) + log(μ/q) - 3/2 ) + 1/q^2 * besselk(0,b₀*q/μ))
# z²Mℂgq₁(q, μ, ζ) = CF/π * CoeffFunc(      z -> (1 + (1-z)^2)/z * b₀/(μ*q) * besselk(1,b₀*q/μ))

function _FUUT_asym(αs::Function, f::Function, D::Function, xB, Q², zh, qT², μ², rtol=_rtol)::Float64
    qT = √qT²
    return αs(μ²)/2π * xB/zh^2 * sum(i -> quark_charge[i]^2 * (
        ( mellinconv(  MCqq₁(qT, √μ², Q²), x -> f(quark_code[i], x, μ²), rtol=rtol)(xB)
        + mellinconv(  MCqg₁(qT, √μ², Q²), x -> f(21,            x, μ²), rtol=rtol)(xB) ) * D(quark_code[i], zh, μ²) +
        ( mellinconv(z²Mℂqq₁(qT, √μ², Q²), z -> D(quark_code[i], z, μ²), rtol=rtol)(zh)
        + mellinconv(z²Mℂgq₁(qT, √μ², Q²), z -> D(21,            z, μ²), rtol=rtol)(zh) ) * f(quark_code[i], xB, μ²)
        ), 1:num_quark)
end

#= (using LHASplitExtend) =========================================================================#

"pid of convolution with splitting kernel"
extpid_1(pid, order)::Int = 1000pid + sign(pid)*( 10 + order )

function FUUT_asym(αs::Function, extf::Function, extD::Function, xB, Q², zh, qT², μ²)::Float64
    xB/(π*zh^2*qT²) * αs(μ²) * sum(i -> quark_charge[i]^2 * (
        1/2π * 2CF * ( - log(qT²/Q²) - 3/2 ) * extf(quark_code[i], xB, μ²) * extD(quark_code[i], zh, μ²) +
        extf(extpid_1(quark_code[i],0), xB, μ²) * extD(quark_code[i], zh, μ²) +
        extf(quark_code[i], xB, μ²) * extD(extpid_1(quark_code[i],0), zh, μ²)
    ), 1:num_quark)
end

#= W+Y construction ===============================================================================#

const _η = 0.34
const _λ = 2/3

Ξ(qT, Q, η=_η, aΞ=8) = exp( - ( qT / ( η*Q ) )^aΞ )
X(qT, λ=_λ, aX=4) = 1 - exp( - ( qT / λ )^aX )

"W+Y construction. `FUUT_tmd` and `FUUT_coll` should follow the interface in `SidisStructFunc`."
function FUUT_WY(FUUT_tmd::Function, FUUT_coll::Function,
        αs::Function, extf::Function, extD::Function,
        xB, Q², zh, qT², μ²,
        η=_η, λ=_λ, rtol=_rtol)::Float64
    qT, Q = √qT², √Q²
    tmd()  = FUUT_tmd(                 xB, Q², zh, qT², μ², rtol)
    coll() = FUUT_coll(                xB, Q², zh, qT², μ², rtol)
    asym() = FUUT_asym(αs, extf, extD, xB, Q², zh, qT², μ²)
    return (
        qT < Q ? Ξ(qT, Q, η) * tmd() : 0
    ) + (
        (qT/Q)^2 > SIDISXSec.Coll.qT²oQ²_low ? X(qT, λ) *( coll() - Ξ(qT, Q, η) * asym() ) : 0
    )
end

"""
    get_sf_WY(data::SidisData, sf_tmd::SidisStructFunc, sf_coll::SidisStructFunc;
        η=_η, λ=_λ)::SidisStructFunc
"""
function get_sf_WY(data::SidisData, sf_tmd::SidisStructFunc, sf_coll::SidisStructFunc;
        η=_η, λ=_λ)::SidisStructFunc
    @info "get_sf_WY: please ensure `data` is the `ZetaSplitExtend` version"
    αs, f, D = data.αs, data.f, data.D
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=_rtol) -> FUUT_WY(sf_tmd.FUUT, sf_coll.FUUT, αs, f, D, xB, Q², zh, qT², μ², η, λ, rtol),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
    )
end

end # module
