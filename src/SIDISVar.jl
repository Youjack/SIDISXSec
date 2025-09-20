#= Part of SIDISXSec.jl                                                                           =#
#= SIDIS variables                                                                                =#

# DIS variables
get_s_from_Ebeam(Ebeam, M=Mp) = 2 * M * Ebeam + M^2
get_s(    xB, y, Q², M=Mp) = Q² /( xB * y ) + M^2
get_xB(s,     y, Q², M=Mp) = Q² /( y * (s - M^2) )
get_y( s, xB,    Q², M=Mp) = Q² /( xB * (s - M^2) )
get_Q²(s, xB, y,     M=Mp) = xB * y * (s - M^2)
get_γ²(   xB,    Q², M=Mp) = (2 * xB * M)^2 / Q²
get_ε(    xB, y, Q², M=Mp) = let γ² = get_γ²(xB, Q², M)
    ( 1 - y - y^2 * γ² /4 )/( 1 - y + y^2 /2 + y^2 * γ² /4 ) end
get_W²(   xB,    Q², M=Mp) = M^2 + (1/xB-1)Q²

# transverse momenta
get_qT²( zh, PhT²) = PhT² / zh^2
get_PhT²(zh, qT² ) = zh^2 * qT²
get_qT²oQ²max(xB, zh) = (1/xB-1) * (1/zh-1)

# spin
get_ST²(d, SL) = d^2 - SL^2

# azimuthal angles
get_ϕ(cosϕ, sinϕ) = angle(cosϕ + im * sinϕ)
get_cos2ϕ(cosϕ, sinϕ) = cosϕ^2 - sinϕ^2
get_sin2ϕ(cosϕ, sinϕ) = 2 * sinϕ * cosϕ
get_trig2ϕ(cosϕ, sinϕ) = ( get_cos2ϕ(cosϕ, sinϕ), get_sin2ϕ(cosϕ, sinϕ) )
get_cos3ϕ(cosϕ, sinϕ) = cosϕ^3 - 3 * cosϕ * sinϕ^2
get_sin3ϕ(cosϕ, sinϕ) = 3 * cosϕ^2 * sinϕ - sinϕ^3
get_trig3ϕ(cosϕ, sinϕ) = ( get_cos3ϕ(cosϕ, sinϕ), get_sin3ϕ(cosϕ, sinϕ) )
get_cosϕ₁₊ϕ₂(cosϕ₁, sinϕ₁, cosϕ₂, sinϕ₂) = cosϕ₁ * cosϕ₂ - sinϕ₁ * sinϕ₂
get_sinϕ₁₊ϕ₂(cosϕ₁, sinϕ₁, cosϕ₂, sinϕ₂) = sinϕ₁ * cosϕ₂ + cosϕ₁ * sinϕ₂
get_cosϕ₁₋ϕ₂(cosϕ₁, sinϕ₁, cosϕ₂, sinϕ₂) = get_cosϕ₁₊ϕ₂(cosϕ₁, sinϕ₁, cosϕ₂, -sinϕ₂)
get_sinϕ₁₋ϕ₂(cosϕ₁, sinϕ₁, cosϕ₂, sinϕ₂) = get_sinϕ₁₊ϕ₂(cosϕ₁, sinϕ₁, cosϕ₂, -sinϕ₂)

"""
Varibles that specify the kinematics of a SIDIS process (on a spin 1/2 target).
- Target mass `M` should never be set to `0`, but `Mh` can be `0`.
- Lepton helicity `λ` should be in `[-1,+1]`. 
- Degree of polarization `d` of the target should be in `[0,1]`, and `SL` should be in `[0,d]`. If
  `ST=0`, or equivalently `|SL|=d`, one can assign `cosϕS=sinϕS=NaN`. For unpolarized target, one
  can assign `d=SL=0` and `cosϕS=sinϕS=NaN`.
- If `PhT²=0`, one can assign `cosϕh=sinϕh=NaN` (`ϕh`-related structure functions vanish due to
  symmetry or oddness when `PhT²→0`). For inclusive DIS, one can assign
  `Mh=zh=cosϕh=sinϕh=PhT²=NaN`.
"""
struct SidisVar
    M     :: Float64
    Mh    :: Float64
    xB    :: Float64
    y     :: Float64
    Q²    :: Float64
    λ     :: Float64
    d     :: Float64
    SL    :: Float64
    cosϕS :: Float64
    sinϕS :: Float64
    zh    :: Float64
    cosϕh :: Float64
    sinϕh :: Float64
    PhT²  :: Float64
end
Base.show(io::IO, var::SidisVar) = let
    M, Mh, xB, y, Q², λ, d, SL, cosϕS, sinϕS, zh, cosϕh, sinϕh, PhT² = exposestruct(var)
    print(io,
        "M = $M GeV, Mh = $Mh GeV\n",
        "√s   = ", √get_s(xB, y, Q², M)       , " GeV \n",
        "xB   = ", xB                         , "     \n",
        "y    = ", y                          , "     \n",
        "Q²   = ", Q²                         , " GeV²\n",
        "λ    = ", λ                          , "     \n",
        "d    = ", d                          , "     \n",
        "SL   = ", SL                         , "     \n",
        "ST   = ", √get_ST²(d, SL)            , "     \n",
        "ϕS   = ", get_ϕ(cosϕS,sinϕS)|>rad2deg, "°    \n",
        "zh   = ", zh                         , "     \n",
        "ϕh   = ", get_ϕ(cosϕh,sinϕh)|>rad2deg, "°    \n",
        "PhT² = ", PhT²                       , " GeV²\n",
    )
end
SidisVar(M, Mh, xB, y, Q², zh, cosϕh, sinϕh, PhT²) =
    SidisVar(M, Mh, xB, y, Q², 0, 0, 0, NaN, NaN, zh, cosϕh, sinϕh, PhT²)
DisVar(M, xB, y, Q², λ, d, SL, cosϕS, sinϕS) =
    SidisVar(M, NaN, xB, y, Q², λ, d, SL, cosϕS, sinϕS, NaN, NaN, NaN, NaN)
DisVar(M, xB, y, Q²) =
    SidisVar(M, NaN, xB, y, Q², 0, 0, 0, NaN, NaN, NaN, NaN, NaN, NaN)
DisVar(v::SidisVar) =
    DisVar(v.M, v.xB, v.y, v.Q², v.λ, v.d, v.SL, v.cosϕS, v.sinϕS)

"""
    above_dis_threshold(var::SidisVar)::Bool

Determine whether `(q+P)² = W² ≥ (M+Mπ)²`.
"""
function above_dis_threshold(var::SidisVar)::Bool
    M, xB, Q² = let v=var; v.M, v.xB, v.Q² end
    return get_W²(xB, Q²) ≥ (M + Mπ)^2
end

"""
    above_sidis_threshold(var::SidisVar)::Bool

Determine whether `(q+P-Ph)² ≥ P²`. `Mh` is set to be at least `Mπ`.
"""
function above_sidis_threshold(var::SidisVar)::Bool
    Mh, xB, Q², zh, PhT² = let v=var; v.Mh, v.xB, v.Q², v.zh, v.PhT² end
    Mh = max(Mh, Mπ)
    γ² = get_γ²(xB, Q²)
    q_dot_Ph = 1/γ² * ( zh * Q² - √( (1 + γ²) * Q² * ( zh^2 * Q² - γ² * (Mh^2 + PhT²) ) ) )
    return 2q_dot_Ph ≤ (1-xB-zh) * Q²/xB + Mh^2
end

include("Frames/lNFrame.jl")

#= QED distortion =================================================================================#

get_ξmin(xB, y, ζ) = (1 - y)/( ζ - xB * y ) # = get_ξmin(xB, y, 0, ζ)
get_ζmin(xB, y   ) = 1 - (1 - xB) * y       # = get_ζmin(xB, y, 0   )
get_ξmin(xB, y, zh, ζ) = max( get_ξmin(xB, y, ζ), y * zh + (1 - y)/ζ   )
get_ζmin(xB, y, zh   ) = max( get_ζmin(xB, y   ), (1 - y)/(1 - y * zh) )

bound0(val) = max(0, val)
bound1(val) = min(val, 1)
boundpm1(val) = clamp(val, -1, +1)

"""
    get_sidis_hat_var(var::SidisVar, ξ, ζ; rot=true)::SidisVar

Get `(ξ,ζ)`-dependent SIDIS variables. `ξ,ζ` are assumed to be in the allowed range.
"""
function get_sidis_hat_var(var::SidisVar, ξ, ζ)::SidisVar
    M, Mh, xB, y, Q², λ, d, SL, cosϕS, sinϕS, zh, cosϕh, sinϕh, PhT² = exposestruct(var)
    ST² = get_ST²(d, SL)
    if !iszero(ST²) && ( isnan(cosϕS) || isnan(sinϕS) )
        throw(ArgumentError("`ϕS=NaN` while `ST≠0` is not well defined."))
    end
    if !iszero(PhT²) && !isnan(PhT²) && ( isnan(cosϕh) || isnan(sinϕh) )
        throw(ArgumentError("`ϕh=NaN` while `PhT²≠0,NaN` is not well defined."))
    end

    # un-hat variables
    γ²  = get_γ²(xB, Q², M)
    lT² = Q² * ( 1 - y - y^2 * γ² /4 )/( y^2 * (1 + γ²) )
    q_dot_S  = 1/γ² * (         - SL * √( (1 + γ²) * Q² * (             γ²                 ) ) )
    q_dot_Ph = 1/γ² * ( zh * Q² -      √( (1 + γ²) * Q² * ( zh^2 * Q² - γ² * (Mh^2 + PhT²) ) ) )
    l_dot_S  = ( iszero(ST² ) ? 0. : - √(lT² * ST² ) * cosϕS ) + ( (1 + y * γ² /2) * q_dot_S                        )/( y * (1 + γ²) )
    l_dot_Ph = ( iszero(PhT²) ? 0. : - √(lT² * PhT²) * cosϕh ) + ( (1 + y * γ² /2) * q_dot_Ph + (1 - y/2) * zh * Q² )/( y * (1 + γ²) )
    ϵ_l_S_P_q  = iszero(ST² ) ? 0. : - √( (1 + 1/γ²) * lT² * ST²  * M^2 * Q² ) * sinϕS
    ϵ_l_Ph_P_q = iszero(PhT²) ? 0. : - √( (1 + 1/γ²) * lT² * PhT² * M^2 * Q² ) * sinϕh
    
    # hat variables
    x̂B = ξ * xB * y /( ξ * ζ - (1 - y) ) |> bound1
    ẑh = ζ * zh * y /( ξ * ζ - (1 - y) ) |> bound1
    ŷ  = 1 - (1 - y)/(ξ * ζ)
    Q̂² = ξ/ζ * Q²
    l̂_dot_S  = ξ * l_dot_S
    l̂_dot_Ph = ξ * l_dot_Ph
    q̂_dot_S  = (ξ - 1/ζ) * l_dot_S  + 1/ζ * q_dot_S
    q̂_dot_Ph = (ξ - 1/ζ) * l_dot_Ph + 1/ζ * q_dot_Ph
    ϵ_l̂_S_P_q̂  = ξ/ζ * ϵ_l_S_P_q
    ϵ_l̂_Ph_P_q̂ = ξ/ζ * ϵ_l_Ph_P_q
    # there should NOT be un-hat variable below !!!
    γ̂²  = get_γ²(x̂B, Q̂², M)
    l̂T² = Q̂² * ( 1 - ŷ - ŷ^2 * γ̂² /4 )/( ŷ^2 * (1 + γ̂²) )
    ŜL  = - q̂_dot_S / √( (1 + 1/γ̂²) * Q̂² ) |> val->clamp(val,-d,d)
    ŜT² = get_ST²(d, ŜL)
    P̂hT² = - Mh^2 + ( - γ̂² * q̂_dot_Ph^2 + 2 * q̂_dot_Ph * ẑh * Q̂² + ẑh^2 * Q̂²^2 )/( (1 + γ̂²) * Q̂² ) |> bound0
    cosϕ̂S = iszero(ST² ) ? NaN : 1/√(l̂T² * ŜT² ) * ( - l̂_dot_S  + ( (1 + ŷ * γ̂² /2) * q̂_dot_S                        )/( ŷ * (1 + γ̂²) ) ) |> boundpm1
    cosϕ̂h = iszero(PhT²) ? NaN : 1/√(l̂T² * P̂hT²) * ( - l̂_dot_Ph + ( (1 + ŷ * γ̂² /2) * q̂_dot_Ph + (1 - ŷ/2) * ẑh * Q̂² )/( ŷ * (1 + γ̂²) ) ) |> boundpm1
    sinϕ̂S = iszero(ST² ) ? NaN : - ϵ_l̂_S_P_q̂  / √( (1 + 1/γ̂²) * l̂T² * ŜT²  * M^2 * Q̂² ) |> boundpm1
    sinϕ̂h = iszero(PhT²) ? NaN : - ϵ_l̂_Ph_P_q̂ / √( (1 + 1/γ̂²) * l̂T² * P̂hT² * M^2 * Q̂² ) |> boundpm1

    return SidisVar(M, Mh, x̂B, ŷ, Q̂², λ, d, ŜL, cosϕ̂S, sinϕ̂S, ẑh, cosϕ̂h, sinϕ̂h, P̂hT²)
end

"The angle from `q` to `q̂`, in degree."
function get_distort_angle(var::SidisVar, ξ, ζ)::Float64
    M, xB, y, Q² = let v=var; v.M, v.xB, v.y, v.Q² end
    
    E  = Q² /( 2 * xB * y * M )
    E′ = (1-y) * E
    
    # 3-vector products
    l²  = E^2
    l′² = E′^2
    ll′ = E * E′ - Q²/2

    # q̂ is always to the left of q, which is proved by -q̂⋅x̂<0 in the Trento frame.
    cosq̂ = ( ξ * l² + 1/ζ * l′² - ( ξ + 1/ζ ) * ll′ )/
        ( √( l² - 2ll′ + l′² ) * √( ξ^2 * l² - 2 * ξ/ζ * ll′ + 1/ζ^2 * l′² ) )
    return rad2deg(acos(cosq̂))
end

#= qed_conv =======================================================================================#

"""
    cumdf(f::Function, If::Function, D::Function, ID::Function,
        var::SidisVar, get_val::Function, bounds::AbstractArray; rtol=_rtol, seed=0)::Vector{Float64}

Get the cumulative distribution function of the quantity returned by `get_val`,
generated by `f(ξ)` and `D(ζ)`.
"""
function cumdf(f::Function, If::Function, D::Function, ID::Function,
        var::SidisVar, get_val::Function, bounds::AbstractArray; rtol=_rtol, seed=0)::Vector{Float64}
    xB, y, zh = let v=var; v.xB, v.y, v.zh end
    ξmin(ζ) = isnan(zh) ? get_ξmin(xB, y, ζ) : get_ξmin(xB, y, zh, ζ)
    ζmin    = isnan(zh) ? get_ζmin(xB, y   ) : get_ζmin(xB, y, zh   )
    probs = Vector{Float64}(undef, length(bounds))
    Threads.@threads for i ∈ eachindex(bounds)
        probs[i] = qed_conv(f, If, D, ID, (ξ,ζ) -> 1,
            (ξ,ζ) -> let v̂al = get_val(get_sidis_hat_var(var, ξ, ζ))
                ζ > ζmin && ξ > ξmin(ζ) && v̂al < bounds[i] + √eps(Float64)
            end, ξmin(1), ζmin, algor=QEDFactorization.INT_MC, rtol=rtol, seed=seed)
    end
    return probs
end
