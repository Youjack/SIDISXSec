#= Part of SIDISXSec.jl                                                                           =#
#= Utilities for calculating SIDIS observables                                                    =#

function change_sidis_var_ϕ(var::SidisVar, ϕS, ϕh)
    ( M, Mh, xB, y, Q², λ, d, SL, _, _, zh, _, _, PhT²,
        γ², ε, lT², ST², qT², q_dot_S, q_dot_Ph,
        _, _, _, _
    ) = exposestruct(var)
    return SidisVar(M, Mh, xB, y, Q², λ, d, SL, cos(ϕS), sin(ϕS), zh, cos(ϕh), sin(ϕh), PhT²,
        γ², ε, lT², ST², qT², q_dot_S, q_dot_Ph)
end

const _trapzNstart = 2
const _trapznmax = 6
const _trapzϕ0 = [
    i * 2π/_trapzNstart
    for i ∈ 0:_trapzNstart-1
]
const _trapzϕ = [
    let N = _trapzNstart * 2^n
        [
            i * 2π/N
            for i ∈ 1:2:N
        ]
    end for n ∈ 1:_trapznmax
]
"Calculate `∫(0,2π) f(ϕ) dϕ` using trapezoidal rule."
function trapzϕ(f, rtol=_rtol)::Float64
    set_prev = mean(f(ϕ) for ϕ ∈ _trapzϕ0)
    set = 0.0
    rerr = 1.0
    for n ∈ 1:_trapznmax
        set_next = mean(f(ϕ) for ϕ ∈ _trapzϕ[n])
        set = ( set_prev + set_next )/2
        if iszero(set)
            @warn "trapzϕ: set = 0"
            return 0.0
        end
        rerr = abs((set - set_prev)/set)
        if rerr ≤ rtol
            return 2π * set
        end
        set_prev = set
    end
    @warn "trapzϕ: rerr = $rerr, sampled points: $(_trapzNstart^(_trapznmax+1))"
    return 2π * set
end
const _trapzϕϕ0 = [
    tuple( i * 2π/_trapzNstart, j * 2π/_trapzNstart )
    for j ∈ 0:_trapzNstart-1, i ∈ 0:_trapzNstart-1
]
const _trapzϕϕ = [
    let N = _trapzNstart * 2^n
        filter(!isnothing, [
            isodd(i) && isodd(j) ? nothing :
            tuple( i * 2π/N , j * 2π/N )
            for j ∈ 0:N-1, i ∈ 0:N-1
        ])
    end for n ∈ 1:_trapznmax
]
"Calculate `∬(0,2π)² f(ϕ₁,ϕ₂) dϕ₁dϕ₂` using trapezoidal rule."
function trapzϕϕ(f, rtol=_rtol)::Float64
    # anisotropic refinement may be used [ChatGPT]
    set_prev = mean(f(ϕ...) for ϕ ∈ _trapzϕϕ0)
    set = 0.0
    rerr = 1.0
    for n ∈ 1:_trapznmax
        set_next = mean(f(ϕ...) for ϕ ∈ _trapzϕϕ[n])
        set = ( set_prev + 3set_next )/4
        rerr = abs((set - set_prev)/set)
        if rerr ≤ rtol
            return (2π)^2 * set
        end
        set_prev = set
    end
    @warn "trapzϕϕ: rerr = $rerr"
    return (2π)^2 * set
end

"""
    SIDIS_mul_xB_Q²_zh_PhT²(data::SidisData, sf::SidisStructFunc,
        var::SidisVar, μ², opt::Options=_opt)::Float64

SIDIS multiplicity `[ dσ /( dxB dQ² dzh dPhT² ) ]/[ dσ /( dxB dQ² ) ]`. \\
Only `FUUL` and `FUUT` are used (without `FLL`).
"""
function SIDIS_mul_xB_Q²_zh_PhT²(data::SidisData, sf::SidisStructFunc,
        var::SidisVar, μ², opt::Options=_opt)::Float64
    var′ = change_sidis_var_ϕ(var, NaN, NaN)
    sf′ = SidisStructFunc(sf.FUUL, sf.FUUT, zerosf, zerosf)
    DIS_xsec = DIS_xsec_xB_Q²_ϕS(data, var′, μ², opt)
    if iszero(DIS_xsec) return 0 end
    SIDIS_xsec = 2π * SIDIS_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf′, var′, μ², opt)
    return SIDIS_xsec / DIS_xsec
end
"""
    SIDISRC_mul_xB_Q²_zh_PhT²(data::SidisData, sf::SidisStructFunc, var::SidisVar, rc::RCData,
        μ², opt::Options=_opt)::Float64

SIDIS multiplicity `[ dσ /( dxB dQ² dzh dPhT² ) ]/[ dσ /( dxB dQ² ) ]` with radiative corrections. \\
(polarizations are not considered)
"""
function SIDISRC_mul_xB_Q²_zh_PhT²(data::SidisData, sf::SidisStructFunc, var::SidisVar, rc::RCData,
        μ², opt::Options=_opt)::Float64
    # DIS
    DIS_xsec = DISRC_xsec_xB_Q²_ϕS(data, var, rc, μ², opt)
    if iszero(DIS_xsec) return 0 end
    # SIDIS
    varϕ(ϕh) = change_sidis_var_ϕ(var, NaN, ϕh)
    if opt.ϕ_algo == ϕGAUSS
        _, fl, Ifl, _, _, Dl, IDl = exposestruct(rc)
        x̂sec(ϕh) = _SIDISRC_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf, varϕ(ϕh), μ², opt, ΣLEPSPIN)
        SIDIS_xsec = 0.0
        SIDIS_xsec +=
            quadgk(ϕh -> 2RCcorner(ξ->Ifl(ξ,μ²), ζ->IDl(ζ,μ²), x̂sec(ϕh), varϕ(ϕh), opt.Mth),
                0,π, rtol=opt.rtol)[1]
        SIDIS_xsec +=
            hcubature(X -> 2RCξedge(ξ->fl(ξ,μ²), ζ->IDl(ζ,μ²), x̂sec(X[1]), varϕ(X[1]), opt.Mth)(X[2]),
                (0,0),(π,1), rtol=opt.rtol, atol=opt.rtol*abs(SIDIS_xsec))[1] +
            hcubature(X -> 2RCζedge(ξ->Ifl(ξ,μ²), ζ->Dl(ζ,μ²), x̂sec(X[1]), varϕ(X[1]), opt.Mth)(X[2]),
                (0,0),(π,1), rtol=opt.rtol, atol=opt.rtol*abs(SIDIS_xsec))[1]
        SIDIS_xsec += opt.incl_rcbulk ?
            hcubature(X -> 2RCbulk(ξ->fl(ξ,μ²), ζ->Dl(ζ,μ²), x̂sec(X[1]), varϕ(X[1]), opt.Mth)(X[2],X[3]),
                (0,0,0),(π,1,1), rtol=opt.rtol, atol=opt.rtol*abs(SIDIS_xsec))[1] : 0.0
    else
        opt′ = Options(opt.rtol/10, opt)
        SIDIS_xsec = trapzϕ(ϕh ->
            SIDISRC_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf, varϕ(ϕh), rc, μ², opt′),
            opt.rtol)
    end    
    return SIDIS_xsec / DIS_xsec
end

"""
    SIDISRC_Aϕh_xB_Q²_zh_PhT²(trig::Function,
        sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt; normalize=true)::Float64

SIDIS azimuthal asymmetry `2⟨trig(ϕh)⟩`. \\
(polarizations are not considered)
"""
function SIDISRC_Aϕh_xB_Q²_zh_PhT²(trig::Function,
        sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt; normalize=true)::Float64
    varϕ(ϕh) = change_sidis_var_ϕ(var, NaN, ϕh)
    opt′ = Options(opt.rtol/10, opt)
    numr = trapzϕ(ϕh ->
        trig(ϕh) * SIDISRC_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf, varϕ(ϕh), rc, μ², opt′),
        opt.rtol)
    denom = normalize ? trapzϕ(ϕh ->
        SIDISRC_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf, varϕ(ϕh), rc, μ², opt′),
        opt.rtol)/2 : 1
    return numr / denom
end
"""
    SIDISRC_AϕSϕh_xB_Q²_zh_PhT²(trig::Function,
        sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt; normalize=true)::Float64

SIDIS azimuthal asymmetry `2⟨trig(ϕS,ϕh)⟩`.
"""
function SIDISRC_AϕSϕh_xB_Q²_zh_PhT²(trig::Function,
        sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt; normalize=true)::Float64
    varϕ(ϕS,ϕh) = change_sidis_var_ϕ(var, ϕS, ϕh)
    opt′ = Options(opt.rtol/10, opt)
    numr = trapzϕ(ϕS -> trapzϕ(ϕh ->
        trig(ϕS,ϕh) * SIDISRC_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf, varϕ(ϕS,ϕh), rc, μ², opt′),
        opt.rtol), opt.rtol)
    denom = normalize ? trapzϕ(ϕS -> trapzϕ(ϕh ->
        SIDISRC_xsec_xB_Q²_ϕS_zh_ϕh_PhT²(sf, varϕ(ϕS,ϕh), rc, μ², opt′),
        opt.rtol), opt.rtol)/2 : 1
    return numr / denom
end
