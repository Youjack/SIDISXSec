module SIDISXSec

include("Constants.jl"); using .Constants
include("QCDData.jl"); using .QCDData
import QuadGK
using QEDFactorization

export exposestruct

export SidisData
export get_toy_data, get_sidis_data, get_sidis_data_from_isospin

export RCData
export get_e⁻_data, get_μ⁻_data

export get_s_from_Ebeam, get_s, get_xB, get_y, get_Q², get_γ², get_ε, get_W²
export get_qT², get_PhT², get_qT²oQ²max
export get_ϕ
export get_ξmin, get_ζmin
export above_dis_threshold, above_sidis_threshold

export SidisVar
export get_sidis_hat_var
export qed_conv

export SidisStructFunc
export get_sf_tmd, get_sf_coll

export DIS_xsec_xB_Q²              , DISRC_xsec_xB_Q²
export SIDIS_xsec_xB_Q²_zh_ϕh_PhT² , SIDISRC_xsec_xB_Q²_zh_ϕh_PhT²
export SIDIS_mul_xB_Q²_zh_PhT²     , SIDISRC_mul_xB_Q²_zh_PhT²
export SIDIS_AUUcosϕh_xB_Q²_zh_PhT², SIDISRC_AUUcosϕh_xB_Q²_zh_PhT²

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
- Collinear PDF `f(id, x, μ²)`
- Collinear FF  `D(id, z, μ²)`
"""
struct SidisData{F1,F2,F3}
    M  :: Float64
    Mh :: Float64
    αs :: F1
    f  :: F2
    D  :: F3
end
Base.show(io::IO, data::SidisData) = print(io,
    "SidisData with M = $(data.M) GeV, Mh = $(data.Mh) GeV"
)

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
    N = A - Z
    return SidisData(
        Float64(M/A), Float64(Mh),
        μ² -> get_αs(fdata, μ²),
        (id, x, μ²) -> ( Z * get_qcd_density(fdata, (-1)^conjN *          id,  μ²)(x)
                       + N * get_qcd_density(fdata, (-1)^conjN * iso_code(id), μ²)(x) )/A,
        (id, z, μ²) ->       get_qcd_density(Ddata, (-1)^conjh *          id,  μ²)(z),
    )
end

#= RCData =========================================================================================#

"""
Data for calculations of radiative corrections.
- `αEM(μ²)`
- `fl(ξ, μ²)`, `Ifl(ξ, μ²)`
- `Dl(ζ, μ²)`, `IDl(ζ, μ²)`
"""
struct RCData{F1,F2,F3,F4,F5}
    αEM :: F1
    fl  :: F2
    Ifl :: F3
    Dl  :: F4
    IDl :: F5
end

function get_e⁻_data()::RCData
    fdata = get_qed_data("QED_electron_DF_nlo")
    Ddata = fdata
    return RCData(
        μ² -> get_αEM(fdata, μ²),
        (ξ, μ²) -> get_qed_density(         fdata, 11, μ²)(ξ),
        (ξ, μ²) -> get_qed_density_integral(fdata, 11, μ²)(ξ),
        (ζ, μ²) -> get_qed_density(         Ddata, 11, μ²)(ζ),
        (ζ, μ²) -> get_qed_density_integral(Ddata, 11, μ²)(ζ),
    )
end
function get_μ⁻_data()::RCData
    fdata = get_qed_data("QED_muon_DF_nlo")
    Ddata = fdata
    return RCData(
        μ² -> get_αEM(fdata, μ²),
        (ξ, μ²) -> get_qed_density(         fdata, 13, μ²)(ξ),
        (ξ, μ²) -> get_qed_density_integral(fdata, 13, μ²)(ξ),
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

"""
    DIS_xsec_xB_Q²(data::SidisData, var::SidisVar, μ², opt::Options=_opt)::Float64

DIS `dσ /( dxB dQ² )/( 2π αEM² )` at LO.
"""
function DIS_xsec_xB_Q²(data::SidisData, var::SidisVar, μ², opt::Options=_opt)::Float64
    M, _, xB, y, Q², _, _, _, _, _, _, _ = exposestruct(var)
    if !above_dis_threshold(var) return 0 end
    if opt.Q_cut^2 > Q²                return 0 end
    if opt.W_cut^2 > get_W²(xB, Q², M) return 0 end
    ε = get_ε(xB, y, Q², M)
    return (y / Q²) *
        2/( y * Q² ) * y^2/2(1-ε) *
        sum(i -> quark_charge[i]^2 * data.f(quark_code[i], xB, μ²), 1:num_quark)
end

"""
    DISRC_xsec_xB_Q²(data::SidisData, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64

DIS `dσ /( dxB dQ² )/( 2π αEM² )` at LO with radiative corrections.
"""
function DISRC_xsec_xB_Q²(data::SidisData, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64
    M, _, xB, y, Q², _, _, _, _, _, _, _ = exposestruct(var)
    x̂sec(ξ,ζ) = let
        v̂ar = get_sidis_hat_var(SidisVar(M, xB, y, Q²), ξ, ζ)
        _, _, _, ŷ, Q̂², _, _, _, _, _, _, _ = exposestruct(v̂ar)
        return (y / Q²)/(ŷ / Q̂²) * 1/ζ * y /( ξ * ζ - (1 - y) ) * # Jacobian
            DIS_xsec_xB_Q²(data, v̂ar, μ², opt)
    end
    _, fl, Ifl, Dl, IDl = exposestruct(rc)
    ξmin(ζ) = get_ξmin(xB, y, ζ)
    ζmin    = get_ζmin(xB, y   )
    #= return qed_conv(
        ξ -> fl(ξ, μ²), ξ -> Ifl(ξ, μ²), ζ -> Dl(ζ, μ²), ζ -> IDl(ζ, μ²), x̂sec,
        (ξ,ζ) -> ζ > ζmin && ξ > ξmin(ζ), ξmin(1), ζmin, rtol=opt.rtol) =#
    fl_x̂sec(ζ) =
        quadgk(ξ -> fl(ξ, μ²) * ( x̂sec(ξ, ζ) - x̂sec(1, ζ) ), ξmin(ζ), 1, rtol=opt.rtol)[1] +
        Ifl(ξmin(ζ), μ²) * x̂sec(1, ζ)
    Dl_fl_x̂sec =
        quadgk(ζ -> Dl(ζ, μ²) * ( fl_x̂sec(ζ) - fl_x̂sec(1) ), ζmin, 1, rtol=opt.rtol)[1] +
        IDl(ζmin, μ²) * fl_x̂sec(1)
    return Dl_fl_x̂sec
end

"""
    SIDIS_xsec_xB_Q²_zh_ϕh_PhT²(sf::SidisStructFunc, var::SidisVar, μ²,
        opt::Options=_opt)::Float64

SIDIS `dσ /( dxB dQ² dzh dϕh dPhT² )/( 2π αEM² )`. \\
(polarizations are not considered)
"""
function SIDIS_xsec_xB_Q²_zh_ϕh_PhT²(sf::SidisStructFunc, var::SidisVar, μ²,
        opt::Options=_opt)::Float64
    M, _, xB, y, Q², _, _, _, zh, cosϕh, sinϕh, PhT² = exposestruct(var)
    if !above_sidis_threshold(var) return 0 end
    if opt.Q_cut^2 > Q²                return 0 end
    if opt.W_cut^2 > get_W²(xB, Q², M) return 0 end
    ε = get_ε(xB, y, Q², M)
    γ² = get_γ²(xB, Q², M)
    qT² = get_qT²(zh, PhT²)
    return (y / Q²) *
        1/( xB * y * Q² ) * y^2/2(1-ε) * (1+γ²/2xB) * (
        ε *                           sf.FUUL(     xB, Q², zh, qT², μ², opt.rtol) +
                                      sf.FUUT(     xB, Q², zh, qT², μ², opt.rtol) +
        ( !isnan(cosϕh) && !isnan(sinϕh) ?
        √(2ε*(1+ε)) * cosϕh *         sf.FUUcosϕh( xB, Q², zh, qT², μ², opt.rtol) +
        ε * get_cos2ϕ(cosϕh, sinϕh) * sf.FUUcos2ϕh(xB, Q², zh, qT², μ², opt.rtol) :
        0 ))
end

"""
    SIDISRC_xsec_xB_Q²_zh_ϕh_PhT²(sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64

SIDIS `dσ /( dxB dQ² dzh dϕh dPhT² )/( 2π αEM² )` with radiative corrections. \\
(polarizations are not considered)
"""
function SIDISRC_xsec_xB_Q²_zh_ϕh_PhT²(sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64
    _, _, xB, y, Q², _, _, _, zh, _, _, _ = exposestruct(var)
    x̂sec(ξ,ζ) = let
        v̂ar = get_sidis_hat_var(var, ξ, ζ)
        _, _, _, ŷ, Q̂², _, _, _, _, _, _, _ = exposestruct(v̂ar)
        return (y / Q²)/(ŷ / Q̂²) * ( y /( ξ * ζ - (1 - y) ) )^2 * # Jacobian
            SIDIS_xsec_xB_Q²_zh_ϕh_PhT²(sf, v̂ar, μ², opt)
    end
    _, fl, Ifl, Dl, IDl = exposestruct(rc)
    ξmin(ζ) = get_ξmin(xB, y, zh, ζ)
    ζmin    = get_ζmin(xB, y, zh   )
    #= return qed_conv(
        ξ -> fl(ξ, μ²), ξ -> Ifl(ξ, μ²), ζ -> Dl(ζ, μ²), ζ -> IDl(ζ, μ²), x̂sec,
        (ξ,ζ) -> ζ > ζmin && ξ > ξmin(ζ), ξmin(1), ζmin, rtol=opt.rtol) =#
    fl_x̂sec(ζ) =
        quadgk(ξ -> fl(ξ, μ²) * ( x̂sec(ξ, ζ) - x̂sec(1, ζ) ), ξmin(ζ), 1, rtol=opt.rtol)[1] +
        Ifl(ξmin(ζ), μ²) * x̂sec(1, ζ)
    Dl_fl_x̂sec =
        quadgk(ζ -> Dl(ζ, μ²) * ( fl_x̂sec(ζ) - fl_x̂sec(1) ), ζmin, 1, rtol=opt.rtol)[1] +
        IDl(ζmin, μ²) * fl_x̂sec(1)
    return Dl_fl_x̂sec
end

#= Observables ====================================================================================#

include("SIDISObservables.jl")

end # module
