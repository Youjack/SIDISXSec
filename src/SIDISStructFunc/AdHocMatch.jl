"""
Part of SIDISStructFunc.jl

An ad hoc matching prescription
"""
module AdHocMatch

using ..SIDISXSec

export get_sf_adhocmatch

squeeze(x) = let k = 10
    ( exp(k*x) - 1 )/( exp(k/2) - 1 )/( 1 + exp(k*(x-1/2)) )
end

function F_adhocmatch(F_tmd::Function, F_coll::Function, xB, Q², zh, qT², μ², rtol=_rtol,
        mode=1)::Float64
    ΔqT² = 0.02^2
    W(qT²) = F_tmd( xB, Q², zh, qT², μ², rtol)
    Z(qT²) = F_coll(xB, Q², zh, qT², μ², rtol)
    if qT²/Q² ≤ SIDISXSec.Coll.qT²oQ²_low + ΔqT²/Q² return max(0.0, W(qT²)) end
    if qT²/Q² ≥ SIDISXSec.TMD.qT²oQ²_up             return max(0.0, Z(qT²)) end
    _Z = Z(qT²); if _Z ≤ 0 return 0.0 end
    _W = W(qT²); if _W ≤ 0 return _Z  end
    if _W > _Z # intermediate region
        if mode == 1
            return _W
        elseif mode == 2
            return _Z
        elseif mode == 3
            return √(_W * _Z)
        elseif mode == 4
            return (_W + _Z) / 2
        elseif mode == 5 # derivative
            _ΔW = _W - W(qT²-ΔqT²)
            _ΔZ = _Z - Z(qT²-ΔqT²)
            λ = _ΔW /( _ΔW + 4_ΔZ ) #|> squeeze
            return _W^(1-λ) * _Z^λ
        elseif mode == 6 # derivative of log
            _ΔW = 1 - W(qT²-ΔqT²) / _W
            _ΔZ = 1 - Z(qT²-ΔqT²) / _Z
            λ = _ΔW /( _ΔW + _ΔZ ) #|> squeeze
            return _W^(1-λ) * _Z^λ
        else
            error("F_adhocmatch: unknown mode = $mode")
        end
    #= else # _W < _Z
        return _ΔW > _ΔZ ? _W : _Z
    end =#
    elseif W(0.0) > 10_Z # large qT region
        return _Z
    else # small qT region
        return _W
    end
end

function get_sf_adhocmatch(sf_tmd::SidisStructFunc, sf_coll::SidisStructFunc;
        mode=1)::SidisStructFunc
    # FUUT only
    return SidisStructFunc(
        (Fargs...) -> F_adhocmatch(sf_tmd.FUUT, sf_coll.FUUT, Fargs..., mode)
    )
end

end # module
