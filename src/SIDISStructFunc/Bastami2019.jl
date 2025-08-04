"""
Part of SIDISStructFunc.jl

[bastami2019semiinclusive] model of unpolarized structure functions
- Following codes in [github.com/prokudin/WW-SIDIS]
- MSTW2008lo + DSS07_pion_plus_lo
"""
module Bastami2019

using ..SIDISXSec
using ..SIDISXSec.QCDData

export get_sf_bastami2019

# unpolarized: first line of Table 1
get_pT²avg() = 0.25
get_PT²avg() = 0.2
get_PhT²avg(zh) = zh^2 * get_pT²avg() + get_PT²avg()

# Boer-Mulders: (Table 4) (A.16) (A.15)
const NBMu = 2.1 * 0.35
const NBMd = -1.111 * -0.9
const αu, βu = 0.73, 3.46
const αd, βd = 1.08, 3.46
const NBMq = (NBMu, NaN, NBMd, NaN)
const αq = (αu, NaN, αd, NaN)
const βq = (βu, NaN, βd, NaN)
𝒩BMq(i, x) = NBMq[i] * x^αq[i]*(1-x)^βq[i] * (αq[i]+βq[i])^(αq[i]+βq[i])/(αq[i]^αq[i]*βq[i]^βq[i])
const MBM = √0.34
get_pT²avgBM() = get_pT²avg() * MBM^2 /( get_pT²avg() + MBM^2 )

# Collins: (Table 3) (A.11) (A.13)
const NCfav = 0.49
const NCdis = -1.00
const NCq = (NCfav, NCdis, NCdis, NCfav) # for π⁺ !!!
const γ = 1.06
const δ = 0.07
𝒩Cq(i, z) = NCq[i] * z^γ*(1-z)^δ * (γ+δ)^(γ+δ)/(γ^γ*δ^δ)
const MC = √1.50
get_PT²avgC() = get_PT²avg() * MC^2 /( get_PT²avg() + MC^2 )

function FUUT(f::Function, D::Function, xB, zh, qT², μ²)::Float64
    PhT² = get_PhT²(zh, qT²)
    PhT²avg = get_PhT²avg(zh)
    return xB * sum(i -> quark_charge[i]^2 *
        f(quark_code[i], xB, μ²) * D(quark_code[i], zh, μ²), 1:num_quark) *
        exp(- PhT² / PhT²avg) /( π * PhT²avg )
end

function FUUcosϕh(f::Function, D::Function, xB, Q², zh, qT², μ²)::Float64
    PhT² = get_PhT²(zh, qT²)
    pT²avg = get_pT²avg()
    PhT²avg = get_PhT²avg(zh)
    # (7.9a)
    cahn = 2xB * sum(i -> quark_charge[i]^2 *
        f(quark_code[i], xB, μ²) * D(quark_code[i], zh, μ²), 1:num_quark) *
        √(qT²/Q²) * (- zh^2 * pT²avg / PhT²avg ) * exp(- PhT² / PhT²avg) /( π * PhT²avg )
    return cahn
end

function FUUcos2ϕh(f::Function, D::Function, xB, zh, qT², μ²)::Float64
    PhT² = get_PhT²(zh, qT²)
    pT²avg = get_pT²avg()
    PT²avg = get_PT²avg()
    pT²avgBM = get_pT²avgBM()
    PhT²avgBM = get_PT²avgC() + zh^2 * get_pT²avgBM()
    # (5.9a) (A.18) (A.12)
    BM = xB * sum(i -> quark_charge[i]^2 *
        𝒩BMq(i, xB) * f(quark_code[i], xB, μ²) *
        𝒩Cq( i, zh) * D(quark_code[i], zh, μ²), (1,3)) * # include only u, d
        -√(ℯ/2) * 1/MBM * pT²avgBM^2/pT²avg * # Boer-Mulders
        √(ℯ/2) * PT²avg * MC^3 /( zh * ( PT²avg + MC^2 )^2 ) * # Collins
        4 * (zh^2*PhT²/PhT²avgBM^2) * exp(-PhT²/PhT²avgBM)/(π*PhT²avgBM)
    return BM
end

function get_sf_bastami2019(data::SidisData)::SidisStructFunc
    f, D = data.f, data.D
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUT(     f, D, xB,     zh, qT², μ²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUcosϕh( f, D, xB, Q², zh, qT², μ²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUcos2ϕh(f, D, xB,     zh, qT², μ²),
    )
end

end # module
