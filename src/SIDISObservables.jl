#= Part of SIDISXsec.jl                                                                           =#
#= SIDIS observables                                                                              =#

get_sidis_var_ϕh(var::SidisVar, ϕh) =
    SidisVar(var.M, var.Mh, var.xB, var.y, var.Q², NaN, NaN, NaN, var.zh,
        cos(ϕh), sin(ϕh), var.PhT²)

"""
    SIDIS_mul_xB_Q²_zh_PhT²(data::SidisData, sf′::SidisStructFunc,
        var::SidisVar, μ², opt::Options=_opt)::Float64

SIDIS multiplicity `[ dσ /( dxB dQ² dzh dPhT² ) ]/[ dσ /( dxB dQ² ) ]`.
Only `FUUL` and `FUUT` are used. \\
(polarizations are not considered)
"""
function SIDIS_mul_xB_Q²_zh_PhT²(data::SidisData, sf′::SidisStructFunc,
        var::SidisVar, μ², opt::Options=_opt)::Float64
    # αEM is not required
    DIS_xsec = DIS_xsec_xB_Q²(data, var, μ², opt)
    if iszero(DIS_xsec) return 0 end
    sf = SidisStructFunc(sf′.FUUL, sf′.FUUT, zerosf, zerosf)
    SIDIS_xsec = 2π * SIDIS_xsec_xB_Q²_zh_ϕh_PhT²(sf, var, μ², opt)
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
    # αEM is not required
    DIS_xsec = DISRC_xsec_xB_Q²(data, var, rc, μ², opt)
    if iszero(DIS_xsec) return 0 end
    SIDIS_xsec = 2 * quadgk( # [0,π] and [-π,0] are symmetric
        ϕh -> SIDISRC_xsec_xB_Q²_zh_ϕh_PhT²(sf, get_sidis_var_ϕh(var, ϕh), rc, μ², opt),
        0, π, rtol=opt.rtol)[1]
    return SIDIS_xsec / DIS_xsec
end

"""
    SIDIS_AUUcosϕh_xB_Q²_zh_PhT²(sf::SidisStructFunc, var::SidisVar, μ²,
        opt::Options=_opt)::Float64

SIDIS asymmetry `AUUcosϕh = 2⟨cosϕh⟩`.
"""
function SIDIS_AUUcosϕh_xB_Q²_zh_PhT²(sf::SidisStructFunc, var::SidisVar, μ²,
        opt::Options=_opt)::Float64
    M, _, xB, y, Q², _, _, _, zh, _, _, PhT² = exposestruct(var)
    ε = get_ε(xB, y, Q², M)
    qT² = get_qT²(zh, PhT²)
    return √(2ε*(1+ε)) * sf.FUUcosϕh( xB, Q², zh, qT², μ², opt.rtol) /
        ( ε * sf.FUUL(xB, Q², zh, qT², μ², opt.rtol)
        +     sf.FUUT(xB, Q², zh, qT², μ², opt.rtol) )
end

"""
    SIDISRC_AUUcosϕh_xB_Q²_zh_PhT²(sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64

SIDIS asymmetry `AUUcosϕh = 2⟨cosϕh⟩` with radiative corrections.
"""
function SIDISRC_AUUcosϕh_xB_Q²_zh_PhT²(sf::SidisStructFunc, var::SidisVar, rc::RCData, μ²,
        opt::Options=_opt)::Float64
    # αEM is not required
    denom = 2 * quadgk( # [0,π] and [-π,0] are symmetric
        ϕh -> SIDISRC_xsec_xB_Q²_zh_ϕh_PhT²(sf, get_sidis_var_ϕh(var, ϕh), rc, μ², opt),
        0, π, rtol=opt.rtol)[1]
    numr = 2 * quadgk( # [0,π] and [-π,0] are symmetric
        ϕh -> cos(ϕh) * SIDISRC_xsec_xB_Q²_zh_ϕh_PhT²(sf, get_sidis_var_ϕh(var, ϕh), rc, μ², opt),
        0, π, rtol=opt.rtol)[1]
    return numr / denom
end
