"""
Part of Coll.jl

Double-parton fragmentation to charged pion. [liu2020power]
"""
module NLPColl

using QuadGK
using HCubature
using ..SIDISXSec
using ..SIDISXSec.Constants, ..SIDISXSec.QCDData

export _get_sf_coll_nlp

const Nc = 3
const Color1 = (Nc^2-1)^2/4Nc^2
const Color8 = (Nc^2-1)/4Nc^3

const fπ = 0.13041
const fK = 0.155

GegenbauerC2(λ, x) = - λ + 2λ * (1+λ) * x^2
ϕπ(x) = 1.81(x*(1-x))^(0.81-0.5) * ( 1 - 0.12GegenbauerC2(0.81, 2x-1) )

get_ŝ(x̂,    Q²        ) = (1-x̂)/x̂ * Q²
get_t̂(   ẑ, Q², qT²oQ²) = - Q² - ẑ * ( -1 + qT²oQ² ) * Q²
get_û(x̂, ẑ, Q²        ) = - ẑ/x̂ * Q²

get_A(ŝ, t̂, û, η) = t̂ + η * (ŝ + û)
get_B(ŝ, t̂, û, η) = ŝ * t̂ - û * (t̂ + η * û)

"(25~28) without `K` factor"
function ℳ²T(C, eq, eq′, ŝ, t̂, û, ξ, ζ)
    ξ̄, ζ̄ = 1-ξ, 1-ζ
    Aξ = get_A(ŝ, t̂, û, ξ)
    Aξ̄ = get_A(ŝ, t̂, û, ξ̄)
    Aζ = get_A(ŝ, t̂, û, ζ)
    Aζ̄ = get_A(ŝ, t̂, û, ζ̄)
    Bξ = get_B(ŝ, t̂, û, ξ)
    Bζ = get_B(ŝ, t̂, û, ζ)
    return (
        # (25)
        + 4C * eq^2 *( - ((t̂+2û)*Aξ̄*Aζ̄+û^2*(t̂-ξ̄*Aζ̄-ζ̄*Aξ̄))/(ξ̄*ζ̄*ŝ^2*Aξ̄*Aζ̄) + (2û^3*(ŝ+t̂+û))/(ŝ*(t̂+û)^2*Aξ̄*Aζ̄) )
        # (26)
        + 4C * eq*eq′ *( ((t̂+2û)*Aξ*Aζ̄+Bξ*Aζ̄-ζ̄*û^2*Aξ)/(ξ*ζ̄*ŝ*û*Aξ*Aζ̄) - (2û*(ŝ+t̂+û)*(t̂+ξ*û))/(ξ*(t̂+û)^2*Aξ*Aζ̄) )
        # (27)
        + 4C * eq*eq′ *( ((t̂+2û)*Aξ̄*Aζ+Bζ*Aξ̄-ξ̄*û^2*Aζ)/(ξ̄*ζ*ŝ*û*Aξ̄*Aζ) - (2û*(ŝ+t̂+û)*(t̂+ζ*û))/(ζ*(t̂+û)^2*Aξ̄*Aζ) )
        # (28)
        + 4C * eq′^2 *( - ((t̂+2û)*Aξ*Aζ+Bξ*Aζ+Bζ*Aξ+ŝ^2*t̂)/(ξ*ζ*û^2*Aξ*Aζ) + (2ŝ*(ŝ+t̂+û)*(t̂+ξ*û)*(t̂+ζ*û))/(ξ*ζ*û*(t̂+û)^2*Aξ*Aζ) )
    )
end
"(32~35) without `K` factor"
function ℳ²L(C, eq, eq′, ŝ, t̂, û, ξ, ζ)
    ξ̄, ζ̄ = 1-ξ, 1-ζ
    Aξ = get_A(ŝ, t̂, û, ξ)
    Aξ̄ = get_A(ŝ, t̂, û, ξ̄)
    Aζ = get_A(ŝ, t̂, û, ζ)
    Aζ̄ = get_A(ŝ, t̂, û, ζ̄)
    return (
        # (32)
        + 16C * eq^2 * (û^3*(ŝ+t̂+û))/(ŝ*(t̂+û)^2*Aξ̄*Aζ̄)
        # (33)
        - 16C * eq*eq′ * (û*(t̂+ξ*û)*(ŝ+t̂+û))/(ξ*(t̂+û)^2*Aξ*Aζ̄)
        # (34)
        - 16C * eq*eq′ * (û*(t̂+ζ*û)*(ŝ+t̂+û))/(ζ*(t̂+û)^2*Aξ̄*Aζ)
        # (35)
        + 16C * eq′^2 * (ŝ*(t̂+ξ*û)*(t̂+ζ*û)*(ŝ+t̂+û))/(ξ*ζ*û*(t̂+û)^2*Aξ*Aζ)
    )
end

const TYPEπ⁺ = 1
const TYPEπ⁻ = 2
const q_code = (
    ( +2, -1 ), # π⁺ with u, d̄
    ( -2, +1 ), # π⁻ with ū, d
)
const q_charge = (
    ( +2/3, +1/3 ), # π⁺ with u, d̄
    ( -2/3, -1/3 ), # π⁻ with ū, d
)
const q′_charge = (
    ( -1/3, -2/3 ), # π⁺ with u, d̄
    ( +1/3, +2/3 ), # π⁻ with ū, d
)

const _rtol = 5e-3
function _C_f_D(αs::Function, f::Function, πtype, ℳ²::Function,
        xB, zh, qT²oQ², Q², rtol)::Float64
    x = xB * ( 1 + zh/(1-zh) * qT²oQ² )
    # if x > 1 return 0 end
    x̂, ẑ = xB/x, zh/1
    ŝ = get_ŝ(x̂,    Q²        )
    t̂ = get_t̂(   ẑ, Q², qT²oQ²)
    û = get_û(x̂, ẑ, Q²        )
    γ² = SIDISXSec.get_γ²(xB, Q², Mp)
    return √(1+γ²) * x̂ * zh/(1-zh) * sum(
        i -> max(1e-10, f(q_code[πtype][i], x, Q²)) *
            fπ^2/16Nc^2 * 16π^2*αs(Q²)^2 *
            # quadgk(ξ -> quadgk(ζ ->
            # ϕπ(ζ) * ϕπ(ξ) * ℳ²(Color1, q_charge[πtype][i], q′_charge[πtype][i], ŝ, t̂, û, ξ, ζ),
            # 0,1, rtol=rtol)[1], 0,1, rtol=rtol)[1],
            hcubature(X -> let ξ = X[1], ζ = X[2]
                ϕπ(ζ) * ϕπ(ξ) * ℳ²(Color1, q_charge[πtype][i], q′_charge[πtype][i], ŝ, t̂, û, ξ, ζ)
            end, (0,0),(1-1e-8,1-1e-8), rtol=rtol)[1],
        (1,2)
    )
end

function _F_coll(αs::Function, f::Function, πtype, ℳ²::Function,
        xB, Q², zh, qT², rtol=_rtol)::Float64
    1/(2π)^3 * xB/(2zh^2*Q²) * _C_f_D(αs, f, πtype, ℳ², xB, zh, qT²/Q², Q², rtol)
end

function _get_sf_coll_nlp(data::SidisData, hcharge)::SidisStructFunc
    αs, f = data.αs, data.f
    πtype::Int = 0
    if     hcharge ≡ +1 πtype = 1
    elseif hcharge ≡ -1 πtype = 2
    else   throw(ArgumentError("Not supported `hcharge`: $hcharge"))
    end
    return SidisStructFunc(
        (xB, Q², zh, qT², μ², rtol=_rtol) -> _F_coll(αs, f, πtype, ℳ²L, xB, Q², zh, qT², rtol),
        (xB, Q², zh, qT², μ², rtol=_rtol) -> _F_coll(αs, f, πtype, ℳ²T, xB, Q², zh, qT², rtol),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
    )
end

end # module
