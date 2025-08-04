"""
Part of SIDISStructFunc.jl

[barone2015phenomenological] model of unpolarized structure functions
- CTEQ6L + DSS07
"""
module Barone2015

using ..SIDISXSec
using ..SIDISXSec.Constants

export get_sf_barone2015

# COMPASS: TABLE I
const A = 0.200
const B = 0.571
const C = 0.031
const MBM = √0.09
const Nd = -1.00
const Nu = -0.45

# HERMES: TABLE II
# const A = 0.126
# const B = 0.506
# const C = 0.037
# const MBM = √0.10
# const Nd = -1.00
# const Nu = -0.49

# Unpolarized: (44~45)
get_pT²avg() = C
get_PT²avg(zh) = A + B * zh^2
get_PhT²avg(zh) = A + (B + C) * zh^2

# Boer-Mulders: (20) (18)
get_pT²avgBM() = get_pT²avg() * MBM^2 /( get_pT²avg() + MBM^2 )
const Nq = (Nu, 0, Nd, 0)
Δf_prefac(i, x) = Nq[i]

# Collins: (46) (24) (22)
const MC = √1.50
get_PT²avgC(zh) = get_PT²avg(zh) * MC^2 /( get_PT²avg(zh) + MC^2 )
const NCfav = 0.49
const NCdis = -1.00
const NCq = (NCfav, NCdis, NCdis, NCfav) # for π⁺ !!!
const γ = 1.06
const δ = 0.07
ΔD_prefac(i, z) = NCq[i] * (γ+δ)^(γ+δ)/(γ^γ*δ^δ) * z^γ*(1-z)^δ

# (32)
get_PhT²avgBM(zh) = get_PT²avgC(zh) + zh^2 * get_pT²avgBM()

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
    pT²avgBM = get_pT²avgBM()
    PT²avg = get_PT²avg(zh)
    PT²avgC = get_PT²avgC(zh)
    PhT²avg = get_PhT²avg(zh)
    PhT²avgBM = get_PhT²avgBM(zh)
    # (26)
    cahn = 2xB * sum(i -> quark_charge[i]^2 *
        f(quark_code[i], xB, μ²) * D(quark_code[i], zh, μ²), 1:4) *
        √(qT²/Q²) * (- zh^2 * pT²avg / PhT²avg ) * exp(- PhT² / PhT²avg) /( π * PhT²avg )
    # (27~28)
    BM = 2xB * ℯ * √(PhT²/Q²) * sum(i -> quark_charge[i]^2 *
        Δf_prefac(i, xB) * f(quark_code[i], xB, μ²) / MBM *
        ΔD_prefac(i, zh) * D(quark_code[i], zh, μ²) / MC, 1:4) *
        exp(- PhT²/PhT²avgBM)/(π*PhT²avgBM^4) *
        pT²avgBM^2*PT²avgC^2/(pT²avg*PT²avg) *
        ( zh^2*pT²avgBM*(PhT²-PhT²avgBM) + PT²avgC*PhT²avgBM )
    return cahn + BM
end

function FUUcos2ϕh(f::Function, D::Function, xB, Q², zh, qT², μ²)::Float64
    PhT² = get_PhT²(zh, qT²)
    pT²avg = get_pT²avg()
    pT²avgBM = get_pT²avgBM()
    PT²avg = get_PT²avg(zh)
    PT²avgC = get_PT²avgC(zh)
    PhT²avg = get_PhT²avg(zh)
    PhT²avgBM = get_PhT²avgBM(zh)
    # (29)
    cahn = 2xB * sum(i -> quark_charge[i]^2 *
        f(quark_code[i], xB, μ²) * D(quark_code[i], zh, μ²), 1:4) *
        qT²/Q² * ( zh^2 * pT²avg / PhT²avg )^2 * exp(- PhT² / PhT²avg) /( π * PhT²avg )
    # (30)
    BM = - xB * ℯ * PhT² * sum(i -> quark_charge[i]^2 *
        Δf_prefac(i, xB) * f(quark_code[i], xB, μ²) / MBM *
        ΔD_prefac(i, zh) * D(quark_code[i], zh, μ²) / MC, 1:4) *
        exp(- PhT²/PhT²avgBM)/(π*PhT²avgBM^3) *
        zh * pT²avgBM^2*PT²avgC^2/(pT²avg*PT²avg)
    return #= cahn + =# BM
end

function get_sf_barone2015(data::SidisData)::SidisStructFunc
    f, D = data.f, data.D
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUT(     f, D, xB,     zh, qT², μ²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUcosϕh( f, D, xB, Q², zh, qT², μ²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUcos2ϕh(f, D, xB, Q², zh, qT², μ²),
    )
end

end # module
