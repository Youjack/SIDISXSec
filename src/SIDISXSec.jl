module SIDISXSec

import QuadGK
using HCubature: hcubature
using Statistics
using QEDFactorization
include("Constants.jl"); using .Constants
include("QCDData.jl"); using .QCDData

export exposestruct

export SidisData
export get_toy_data, get_sidis_data, get_polsidis_data

export RCData
export get_e⁻_data, get_μ⁻_data

export get_s_from_Ebeam, get_s, get_xB, get_y, get_Q², get_γ², get_ε, get_W²
export get_qT², get_PhT², get_qT²oQ²max
export get_ϕ
export get_ξmin, get_ζmin
export above_dis_threshold, above_sidis_threshold

export SidisVar, DisVar
export get_sidis_hat_var
export qed_conv

export SidisStructFunc
export get_sf_tmd, get_sf_coll

export DIS_xsec_xB_Q²_ϕS             , DISRC_xsec_xB_Q²_ϕS
export SIDIS_xsec_xB_Q²_ϕS_zh_ϕh_PhT², SIDISRC_xsec_xB_Q²_ϕS_zh_ϕh_PhT²
export SIDIS_mul_xB_Q²_zh_PhT²       , SIDISRC_mul_xB_Q²_zh_PhT²
export trapzϕ, trapzϕϕ
export SIDISRC_Aϕh_xB_Q²_zh_PhT², SIDISRC_AϕSϕh_xB_Q²_zh_PhT²

const _rtol = 5e-3
const _order = 7

@generated function exposestruct(S)
    fields = fieldnames(S)
    return Expr(:tuple, [:(S.$field) for field ∈ fields]...)
end

function quadgk(f, a, b; rtol=_rtol)::Tuple{Float64,Float64}
    # `quadgk` also calls `handle_infinities` !!!
    QuadGK.do_quadgk(f, (Float64(a),Float64(b)),
        7, nothing, rtol, 10^7, QuadGK.norm, nothing, nothing)
end

#= SidisData ======================================================================================#

"""
Non-perturbative data for SIDIS calculations
- `M` (target mass), `Mh`
- `αs(μ²)`
- Collinear PDF          `f(id, x, μ²)`
- Collinear helicity PDF `g(id, x, μ²)`
- Collinear FF           `D(id, z, μ²)`
"""
struct SidisData{F1,F2,F3,F4}
    M  :: Float64
    Mh :: Float64
    αs :: F1
    f  :: F2
    g  :: F3
    D  :: F4
end
Base.show(io::IO, data::SidisData) = print(io,
    "SidisData with M = $(data.M) GeV, Mh = $(data.Mh) GeV"
)
zeropdf(id, x, μ²) = 0.0
SidisData(M, Mh, αs, f, D) = SidisData(M, Mh, αs, f, zeropdf, D)

function get_toy_data()::SidisData
    return SidisData(
        Mp, 0.0,
        μ² -> 0.118,
        (id, x, μ²) -> get_pdf_toy(id)(x),
        (id, z, μ²) -> get_pdf_toy(id)(z),
    )
end
function get_sidis_data(M, fname::String, Mh, Dname::String;
        A=1, Z=1, conjN=false, conjh=false)::SidisData
    fdata = get_qcd_data(fname)
    Ddata = get_qcd_data(Dname)
    return SidisData(
        Float64(M/A), Float64(Mh),
        μ² -> get_αs(fdata, μ²),
        (id, x, μ²) -> get_qcd_density(fdata, id, μ², conjN, A, Z)(x),
        (id, z, μ²) -> get_qcd_density(Ddata, id, μ², conjh)(z),
    )
end
function get_polsidis_data(M, fname::String, gname::String, Mh, Dname::String;
        A=1, Z=1, conjN=false, conjh=false)::SidisData
    fdata = get_qcd_data(fname)
    gdata = get_qcd_data(gname)
    Ddata = get_qcd_data(Dname)
    return SidisData(
        Float64(M/A), Float64(Mh),
        μ² -> get_αs(fdata, μ²),
        (id, x, μ²) -> get_qcd_density(fdata, id, μ², conjN, A, Z)(x),
        (id, x, μ²) -> get_qcd_density(gdata, id, μ², conjN, A, Z)(x),
        (id, z, μ²) -> get_qcd_density(Ddata, id, μ², conjh)(z),
    )
end

#= RCData =========================================================================================#

"""
Data for calculations of radiative corrections.
- `αEM(μ²)`
- `fl(ξ, μ²)`, `Ifl(ξ, μ²)`
- `gl(ξ, μ²)`, `Igl(ξ, μ²)`
- `Dl(ζ, μ²)`, `IDl(ζ, μ²)`
"""
struct RCData{F1,F2,F3,F4,F5,F6,F7}
    αEM :: F1
    fl  :: F2
    Ifl :: F3
    gl  :: F4
    Igl :: F5
    Dl  :: F6
    IDl :: F7
end
const zerorc = RCData(
    μ² -> 0.0,
    (ξ, μ²) -> 0.0,
    (ξ, μ²) -> 1.0,
    (ξ, μ²) -> 0.0,
    (ξ, μ²) -> 1.0,
    (ζ, μ²) -> 0.0,
    (ζ, μ²) -> 1.0,
)

function get_e⁻_data()::RCData
    fdata = get_qed_data("QED_electron_DF_nlo")
    gdata = get_qed_data("QED_electron_helicityDF_nlo")
    Ddata = fdata
    return RCData(
        μ² -> get_αEM(fdata, μ²),
        (ξ, μ²) -> get_qed_density(         fdata, 11, μ²)(ξ),
        (ξ, μ²) -> get_qed_density_integral(fdata, 11, μ²)(ξ),
        (ξ, μ²) -> get_qed_density(         gdata, 11, μ²)(ξ),
        (ξ, μ²) -> get_qed_density_integral(gdata, 11, μ²)(ξ),
        (ζ, μ²) -> get_qed_density(         Ddata, 11, μ²)(ζ),
        (ζ, μ²) -> get_qed_density_integral(Ddata, 11, μ²)(ζ),
    )
end
function get_μ⁻_data()::RCData
    fdata = get_qed_data("QED_muon_DF_nlo")
    gdata = get_qed_data("QED_muon_helicityDF_nlo")
    Ddata = fdata
    return RCData(
        μ² -> get_αEM(fdata, μ²),
        (ξ, μ²) -> get_qed_density(         fdata, 13, μ²)(ξ),
        (ξ, μ²) -> get_qed_density_integral(fdata, 13, μ²)(ξ),
        (ξ, μ²) -> get_qed_density(         gdata, 13, μ²)(ξ),
        (ξ, μ²) -> get_qed_density_integral(gdata, 13, μ²)(ξ),
        (ζ, μ²) -> get_qed_density(         Ddata, 13, μ²)(ζ),
        (ζ, μ²) -> get_qed_density_integral(Ddata, 13, μ²)(ζ),
    )
end

#= SIDIS variables ================================================================================#

include("SIDISVar.jl")

#= SIDIS structure functions ======================================================================#

include("SIDISStructFunc.jl")

#= Cross-sections =================================================================================#

"Options for evaluating SIDIS xsec"
struct Options
    rtol  :: Float64
    Q_cut :: Float64
    W_cut :: Float64
end
Options(;
    rtol  = _rtol,
    Q_cut = 0.0,
    W_cut = 0.0
) = Options(rtol, Q_cut, W_cut)
const _opt = Options()

const ΣLEPSPIN = 1
const ΔLEPSPIN = 2

function _DIS_xsec_xB_Q²_ϕS(data::SidisData, var::SidisVar, μ², opt::Options=_opt,
        lepspin_mode::Int=0)::Float64
    M, xB, y, Q², λ, SL = let v=var; v.M, v.xB, v.y, v.Q², v.λ, v.SL end
    if !above_dis_threshold(var) return 0 end
    if opt.Q_cut^2 > Q²                return 0 end
    if opt.W_cut^2 > get_W²(xB, Q², M) return 0 end
    ε = get_ε(xB, y, Q², M)
    return (y / Q²) *
        2/( y * Q² ) * y^2/2(1-ε) *
        sum(i -> quark_charge[i]^2 * (
            + ( lepspin_mode == ΔLEPSPIN ? 0 :
                data.f(quark_code[i], xB, μ²)
            )
            + ( lepspin_mode == ΣLEPSPIN || iszero(λ) || iszero(SL) ? 0 :
                (lepspin_mode == ΔLEPSPIN ? 1 : λ) *
                SL * √(1-ε^2) * data.g(quark_code[i], xB, μ²)
            )
        ), 1:num_quark)
end
"""
    DIS_xsec_xB_Q²_ϕS(data::SidisData, var::SidisVar, μ², opt::Options=_opt)::Float64

DIS `dσ /( dxB dQ² dϕS )/ αEM²` at LO. (`ϕS`-independent at leading twist)
"""
function DIS_xsec_xB_Q²_ϕS(data::SidisData, var::SidisVar, μ², opt::Options=_opt)::Float64
    return _DIS_xsec_xB_Q²_ϕS(data, var, μ², opt)
end

function _SIDIS_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf::SidisStructFunc, var::SidisVar, μ²,
        opt::Options=_opt, lepspin_mode::Int=0)::Float64
    M, _, xB, y, Q², λ, d, SL, cosϕS, sinϕS, zh, cosϕh, sinϕh, PhT² = exposestruct(var)
    if !above_sidis_threshold(var) return 0 end
    if opt.Q_cut^2 > Q²                return 0 end
    if opt.W_cut^2 > get_W²(xB, Q², M) return 0 end
    qT² = get_qT²(zh, PhT²)
    Fargs = (xB, Q², zh, qT², μ², opt.rtol)
    ε = get_ε(xB, y, Q², M)
    γ² = get_γ²(xB, Q², M)
    ST = √get_ST²(d, SL)
    return (y / Q²) * 1/( xB * y * Q² ) * y^2/2(1-ε) * (1+γ²/2xB) * (
        + ( lepspin_mode == ΔLEPSPIN ? 0 :
        + (
            + ε * sf.FUUL(Fargs...)
            +     sf.FUUT(Fargs...)
        )
        + ( isnan(cosϕh) || isnan(sinϕh) ? 0 :
            + √(2ε*(1+ε)) * cosϕh                  * sf.FUUcosϕh( Fargs...)
            + ε           * get_cos2ϕ(cosϕh,sinϕh) * sf.FUUcos2ϕh(Fargs...)
        )
        # + ( iszero(SL) || isnan(cosϕh) || isnan(sinϕh) ? 0 : SL * (
        #     + ε           * get_sin2ϕ(cosϕh,sinϕh) * sf.FULsin2ϕh(Fargs...)
        #     + √(2ε*(1+ε)) * sinϕh                  * sf.FULsinϕh( Fargs...)
        # ))
        + ( iszero(ST) || isnan(cosϕS) || isnan(sinϕS) ? 0 : ST * (
            # + √(2ε*(1+ε)) * sinϕS * sf.FUTsinϕS(Fargs...)
        + ( isnan(cosϕh) || isnan(sinϕh) ? 0 :
            # + ε           * get_sinϕ₁₋ϕ₂(cosϕh,sinϕh,cosϕS,sinϕS)                * sf.FUTLsinϕh₋ϕS(Fargs...)
            +               get_sinϕ₁₋ϕ₂(cosϕh,sinϕh,cosϕS,sinϕS)                * sf.FUTTsinϕh₋ϕS(Fargs...)
            + ε           * get_sinϕ₁₊ϕ₂(cosϕh,sinϕh,cosϕS,sinϕS)                * sf.FUTsinϕh₊ϕS( Fargs...)
            # + √(2ε*(1+ε)) * get_sinϕ₁₋ϕ₂(get_trig2ϕ(cosϕh,sinϕh)...,cosϕS,sinϕS) * sf.FUTsin2ϕh₋ϕS(Fargs...)
            + ε           * get_sinϕ₁₋ϕ₂(get_trig3ϕ(cosϕh,sinϕh)...,cosϕS,sinϕS) * sf.FUTsin3ϕh₋ϕS(Fargs...)
        )
        ))
        )
        + ( lepspin_mode == ΣLEPSPIN || iszero(λ) ? 0 : (lepspin_mode == ΔLEPSPIN ? 1 : λ) * (
        # + ( isnan(cosϕh) || isnan(sinϕh) ? 0 :
        #     + √(2ε*(1-ε)) * sinϕh * FLUsinϕh(Fargs...)
        # )
        + ( iszero(SL) ? 0 : SL * (
            + √(1-ε^2) * sf.FLL(Fargs...)
        # + ( isnan(cosϕh) || isnan(sinϕh) ? 0 :
        #     + √(2ε*(1-ε)) * cosϕh * FLLcosϕh(Fargs...)
        # )
        ))
        # + ( iszero(ST) || isnan(cosϕS) || isnan(sinϕS) ? 0 :
        #     + √(2ε*(1-ε)) * cosϕS * FLTcosϕS(Fargs...)
        # + ( isnan(cosϕh) || isnan(sinϕh) ? 0 :
        #     + √(1-ε^2)    * get_cosϕ₁₋ϕ₂(cosϕh,sinϕh,cosϕS,sinϕS)                * FLTcosϕh₋ϕS( Fargs...)
        #     + √(2ε*(1-ε)) * get_cosϕ₁₋ϕ₂(get_trig2ϕ(cosϕh,sinϕh)...,cosϕS,sinϕS) * FLTcos2ϕh₋ϕS(Fargs...)
        # )
        # )
        ))
    )
end
"""
    SIDIS_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf::SidisStructFunc, var::SidisVar, μ²,
        opt::Options=_opt)::Float64

SIDIS `dσ /( dxB dQ² dϕS dzh dϕh dPhT² )/ αEM²`.
- One can assign `ϕS=NaN` or `ϕh=NaN` if the corresponding asymmetry is not needed.
"""
function SIDIS_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf::SidisStructFunc, var::SidisVar, μ²,
        opt::Options=_opt)::Float64
    return _SIDIS_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf, var, μ², opt)
end

#= QED-radiation corrected cross-sections =========================================================#

function qed_conv_xsec(fl, Ifl, Dl, IDl, x̂sec, y, x, z, rtol)
    ξmx = (1-y)/(1-x*y)
    ζmz = (1-y)/(1-z*y)
    ζmx = 1-y+x*y
    ξmz = 1-y+z*y
    ξmxz = max(ξmx, ξmz)
    ζmxz = max(ζmx, ζmz)
    jac_x̂sec′(ξ̃m, ζ̃m) = let
        A = (ξ̃m-ζmz)/(1-ζmz)
        B = (ζ̃m-ξmx)/(1-ξmx)
        C = (1-ξ̃m)/(1-ζmz)*ζmz
        D = (1-ζ̃m)/(1-ξmx)*ξmx
        Δ = √(A^2*B^2+(C-D)^2+2A*B*(C+D))
        ξ = (A*B+C-D+Δ)/2B |> bound1
        ζ = (A*B-C+D+Δ)/2A |> bound1
        E = (ξ-ξmx)/(1-ξmx)
        F = (ζ-ζmz)/(1-ζmz)
        jac = ξ*ζ*E*F/(ξ^2*ζ^2-C*D)
        return jac * fl(ξ) * Dl(ζ) * ( x̂sec(ξ,ζ) - x̂sec(ξ,1.) - x̂sec(1.,ζ) + x̂sec(1.,1.) )
    end
    return (
        + hcubature(X -> jac_x̂sec′(X[1],X[2]), (ξmxz,ζmxz),(1.,1.), rtol=rtol)[1]
        + quadgk(ξ -> fl(ξ) * ( x̂sec(ξ,1.) - x̂sec(1.,1.) ), ξmxz,1., rtol=rtol)[1] * IDl(ζmxz)
        + quadgk(ζ -> Dl(ζ) * ( x̂sec(1.,ζ) - x̂sec(1.,1.) ), ζmxz,1., rtol=rtol)[1] * Ifl(ξmxz)
        + Ifl(ξmxz) * IDl(ζmxz) * x̂sec(1.,1.)
    )
end
#= function qed_conv_xsec(fl, Ifl, Dl, IDl, x̂sec, y, x, z, rtol)
    ξmin(ζ) = get_ξmin(x, y, z, ζ)
    ζmin    = get_ζmin(x, y, z   )
    return qed_conv(fl, Ifl, Dl, IDl, x̂sec,
        (ξ,ζ) -> ζ > ζmin && ξ > ξmin(ζ), ξmin(1), ζmin, rtol=rtol)
end =#
#= function qed_conv_xsec(fl, Ifl, Dl, IDl, x̂sec, y, x, z, rtol)
    ξmin(ζ) = get_ξmin(x, y, z, ζ)
    ζmin    = get_ζmin(x, y, z   )
    fl_x̂sec(ζ) =
        quadgk(ξ -> fl(ξ) * ( x̂sec(ξ, ζ) - x̂sec(1, ζ) ), ξmin(ζ), 1, rtol=rtol)[1] +
        Ifl(ξmin(ζ)) * x̂sec(1, ζ)
    Dl_fl_x̂sec =
        quadgk(ζ -> Dl(ζ) * ( fl_x̂sec(ζ) - fl_x̂sec(1) ), ζmin, 1, rtol=rtol)[1] +
        IDl(ζmin) * fl_x̂sec(1)
    return Dl_fl_x̂sec
end =#

"""
    DISRC_xsec_xB_Q²_ϕS(data::SidisData, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64

DIS `dσ /( dxB dQ² dϕS )/ αEM²` at LO with radiative corrections.
"""
function DISRC_xsec_xB_Q²_ϕS(data::SidisData, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64
    xB, y, Q², λ = let v=var; v.xB, v.y, v.Q², v.λ end
    x̂sec(ξ,ζ, lepspin_mode) = let
        v̂ar = get_sidis_hat_var(DisVar(var), ξ, ζ)
        ŷ, Q̂² = let v=v̂ar; v.y, v.Q² end
        return (y / Q²)/(ŷ / Q̂²) * 1/ζ * y /( ξ * ζ - (1 - y) ) * # Jacobian
            _DIS_xsec_xB_Q²_ϕS(data, v̂ar, μ², opt, lepspin_mode)
    end
    Σx̂sec(ξ,ζ) = x̂sec(ξ,ζ, ΣLEPSPIN)
    Δx̂sec(ξ,ζ) = x̂sec(ξ,ζ, ΔLEPSPIN)
    _, fl, Ifl, gl, Igl, Dl, IDl = exposestruct(rc)
    return (
        + qed_conv_xsec(ξ->fl(ξ,μ²), ξ->Ifl(ξ,μ²), ζ->Dl(ζ,μ²), ζ->IDl(ζ,μ²),
            Σx̂sec, y, xB, 0., opt.rtol)
        + ( iszero(λ) ? 0 : λ*
          qed_conv_xsec(ξ->gl(ξ,μ²), ξ->Igl(ξ,μ²), ζ->Dl(ζ,μ²), ζ->IDl(ζ,μ²),
            Δx̂sec, y, xB, 0., opt.rtol) )
    )
end

"""
    SIDISRC_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64

SIDIS `dσ /( dxB dQ² dϕS dzh dϕh dPhT² )/ αEM²` with radiative corrections.
"""
function SIDISRC_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64
    xB, y, Q², λ, zh = let v=var; v.xB, v.y, v.Q², v.λ, v.zh end
    x̂sec(ξ,ζ, lepspin_mode) = let
        v̂ar = get_sidis_hat_var(var, ξ, ζ)
        ŷ, Q̂² = let v=v̂ar; v.y, v.Q² end
        return (y / Q²)/(ŷ / Q̂²) * ( y /( ξ * ζ - (1 - y) ) )^2 * # Jacobian
            _SIDIS_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf, v̂ar, μ², opt, lepspin_mode)
    end
    Σx̂sec(ξ,ζ) = x̂sec(ξ,ζ, ΣLEPSPIN)
    Δx̂sec(ξ,ζ) = x̂sec(ξ,ζ, ΔLEPSPIN)
    _, fl, Ifl, gl, Igl, Dl, IDl = exposestruct(rc)
    return (
        + qed_conv_xsec(ξ->fl(ξ,μ²), ξ->Ifl(ξ,μ²), ζ->Dl(ζ,μ²), ζ->IDl(ζ,μ²),
            Σx̂sec, y, xB, zh, opt.rtol)
        + ( iszero(λ) ? 0 : λ*
          qed_conv_xsec(ξ->gl(ξ,μ²), ξ->Igl(ξ,μ²), ζ->Dl(ζ,μ²), ζ->IDl(ζ,μ²),
            Δx̂sec, y, xB, zh, opt.rtol) )
    )
end

#= Observables ====================================================================================#

include("SIDISObservables.jl")

end # module
