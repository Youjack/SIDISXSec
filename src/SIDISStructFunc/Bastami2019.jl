"""
Part of SIDISStructFunc.jl

[bastami2019semiinclusive] model of unpolarized structure functions
- Following codes in [github.com/prokudin/WW-SIDIS]
- MSTW2008lo + DSS07_pion_plus_lo
- Support p/n target with π⁺/π⁻ in final states
"""
module Bastami2019

using ..SIDISXSec
using ..SIDISXSec.QCDData

export get_sf_bastami2019

#= TMDs ===========================================================================================#

# masses (actually unused)
const MN = 1
const Mh = 1

# Gaussian ansatz
get_PhT²avg(zh, pT²avg, PT²avg) = zh^2 * pT²avg + PT²avg
Gauss(PhT², PhT²avg) = exp( - PhT² / PhT²avg )/( π * PhT²avg )

# Unpolarized: first line of Table 1
const pT²avg = 0.25
const PT²avg = 0.2
_fperp1(x) = 1/x * pT²avg/2MN^2

# Helicity: section A.2
const pT²avgg = 0.76pT²avg

# Transversity: (Table 3) (A.7) (A.9)
const NTu = 0.46
const NTd = -1.0
const NT = (NTu, NaN, NTd, NaN)
const αT = 1.11
const βT = 3.64
𝒩T(i, x) = NT[i] * x^αT*(1-x)^βT * (αT+βT)^(αT+βT)/(αT^αT*βT^βT)
const pT²avgT = 0.25

# Boer-Mulders: (Table 4) (A.16) (A.15) (A.18)
const NBMu = 2.1 * 0.35
const NBMd = -1.111 * -0.9
const αBMu, βBMu = 0.73, 3.46
const αBMd, βBMd = 1.08, 3.46
const NBM = (NBMu, NaN, NBMd, NaN)
const αBM = (αBMu, NaN, αBMd, NaN)
const βBM = (βBMu, NaN, βBMd, NaN)
𝒩BM(i, x) = NBM[i] * x^αBM[i]*(1-x)^βBM[i] * (αBM[i]+βBM[i])^(αBM[i]+βBM[i])/(αBM[i]^αBM[i]*βBM[i]^βBM[i])
const MBM = √0.34
const pT²avgBM = pT²avg * MBM^2 /( pT²avg + MBM^2 )
_hperp1(i, x) = -√(ℯ/2) * 1/(MN*MBM) * pT²avgBM^2/pT²avg * 𝒩BM(i, x)

# Sivers: (Table 2) (A.4) (A.3) (A.1)
const NSu = 0.4
const NSd = -0.97
const NS = (NSu, NaN, NSd, NaN)
const αSu, αSd = 0.35, 0.44
const βSu, βSd = 2.6, 0.9
const αS = (αSu, NaN, αSd, NaN)
const βS = (βSu, NaN, βSd, NaN)
𝒩S(i, x) = NS[i] * x^αS[i]*(1-x)^βS[i] * (αS[i]+βS[i])^(αS[i]+βS[i])/(αS[i]^αS[i]*βS[i]^βS[i])
const M1 = √0.19
const pT²avgS = pT²avg * M1^2 /( pT²avg + M1^2 )
_fTperp1(i, x) = -√(ℯ/2) * 1/(MN*M1) * pT²avgS^2/pT²avg * 𝒩S(i, x)

# Collins: (Table 3) (A.11) (A.12) (A.13)
const NCfav = 0.49
const NCdis = -1.00
const NC = (
    (NCfav, NCdis, NCdis, NCfav), # π⁺
    (NCdis, NCfav, NCfav, NCdis), # π⁻
)
const γC = 1.06
const δC = 0.07
𝒩C(type, i, z) = NC[type][i] * z^γC*(1-z)^δC * (γC+δC)^(γC+δC)/(γC^γC*δC^δC)
const MC = √1.50
const PT²avgC = PT²avg * MC^2 /( PT²avg + MC^2 )
_Hperp1(type, i, z) = √(ℯ/2) * PT²avg * MC^3 /( z * Mh * ( MC^2 + PT²avg )^2 ) * 𝒩C(type, i, z)

#= Structure Functions ============================================================================#

get_πtype(πcharge) = πcharge == +1 ? 1 : ( πcharge == -1 ? 2 :
    throw(ArgumentError("Unknow πcharge = $πcharge")) )

function FUUT(f, D, xB, zh, qT², μ²)::Float64
    PhT² = get_PhT²(zh, qT²)
    PhT²avg = get_PhT²avg(zh, pT²avg, PT²avg)
    return xB * sum(i -> quark_charge[i]^2
        * f(quark_code[i], xB, μ²)
        * D(quark_code[i], zh, μ²)
        * Gauss(PhT²,PhT²avg),
    1:num_quark)
end

function FUUcosϕh(f, D, xB, Q², zh, qT², μ²)::Float64
    PhT² = get_PhT²(zh, qT²)
    PhT²avg = get_PhT²avg(zh, pT²avg, PT²avg)
    # (7.9a)
    cahn = 2MN/√Q² * xB * sum(i -> quark_charge[i]^2
        * -xB*_fperp1(xB) * f(quark_code[i], xB, μ²)
        * D(quark_code[i], zh, μ²)
        * 2MN * (zh*√PhT²/PhT²avg) * Gauss(PhT²,PhT²avg),
    1:num_quark)
    return cahn
end

"Boer-Mulders asymmetry"
function FUUcos2ϕh(f, D, πcharge, xB, zh, qT², μ²)::Float64
    type = get_πtype(πcharge)
    PhT² = get_PhT²(zh, qT²)
    PhT²avg = get_PhT²avg(zh, pT²avgBM, PT²avgC)
    # (5.9a)
    BM = xB * sum(i -> quark_charge[i]^2
        * _hperp1(i, xB) * f(quark_code[i], xB, μ²)
        * _Hperp1(type, i, zh) * D(quark_code[i], zh, μ²)
        * 4MN*Mh * (zh^2*PhT²/PhT²avg^2) * Gauss(PhT²,PhT²avg),
    (1,3)) # include only u, d
    return BM
end

"Sivers asymmetry"
function FUTTsinϕh₋ϕS(f, D, xB, zh, qT², μ²)::Float64
    PhT² = get_PhT²(zh, qT²)
    PhT²avg = get_PhT²avg(zh, pT²avgS, PT²avg)
    # (5.7a)
    sivers = - xB * sum(i -> quark_charge[i]^2
        * _fTperp1(i, xB) * f(quark_code[i], xB, μ²)
        * D(quark_code[i], zh, μ²)
        * 2MN * (zh*√PhT²/PhT²avg) * Gauss(PhT²,PhT²avg),
    (1,3)) # include only u, d
    return sivers
end

"Collins asymmetry"
function FUTsinϕh₊ϕS(f, g, D, πcharge, xB, zh, qT², μ²)::Float64
    type = get_πtype(πcharge)
    PhT² = get_PhT²(zh, qT²)
    PhT²avg = get_PhT²avg(zh, pT²avgT, PT²avgC)
    # (5.8a)
    collins = xB * sum(i -> quark_charge[i]^2
        * 1/2 * 𝒩T(i, xB) * ( f(quark_code[i], xB, μ²) + g(quark_code[i], xB, μ²) ) # (A.7)
        * _Hperp1(type, i, zh) * D(quark_code[i], zh, μ²)
        * 2Mh * (zh*√PhT²/PhT²avg) * Gauss(PhT²,PhT²avg),
    (1,3)) # include only u, d
    return collins
end

function FLL(g, D, xB, zh, qT², μ²)::Float64
    PhT² = get_PhT²(zh, qT²)
    PhT²avg = get_PhT²avg(zh, pT²avgg, PT²avg)
    return xB * sum(i -> quark_charge[i]^2
        * g(quark_code[i], xB, μ²)
        * D(quark_code[i], zh, μ²)
        * Gauss(PhT²,PhT²avg),
    1:num_quark)
end

function get_sf_bastami2019(data::SidisData; πcharge=+1)::SidisStructFunc
    f, g, D = data.f, data.g, data.D
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUT(        f,    D,          xB,     zh, qT², μ²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUcosϕh(    f,    D,          xB, Q², zh, qT², μ²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUcos2ϕh(   f,    D, πcharge, xB,     zh, qT², μ²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUTTsinϕh₋ϕS(f,    D,          xB,     zh, qT², μ²),
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUTsinϕh₊ϕS( f, g, D, πcharge, xB,     zh, qT², μ²),
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=0.0) -> FLL(            g, D,          xB,     zh, qT², μ²),
    )
end

end # module
