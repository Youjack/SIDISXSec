"""
Part of SIDISStructFunc.jl

[anselmino2014unpolarised] model of unpolarized structure functions
"""
module Anselmino2014

using ..SIDISXSec
using ..SIDISXSec.QCDData

export get_sf_anselmino2014

# (Table 1) HERMES (z < 0.6)
const pT²avg = 0.57
const PT²avg = 0.12

function FUUT(f::Function, D::Function, xB, zh, qT², μ²)::Float64
    PhT² = get_PhT²(zh, qT²)
    PhT²avg = zh^2 * pT²avg + PT²avg
    return xB * sum(i -> quark_charge[i]^2 *
        f(quark_code[i], xB, μ²) * D(quark_code[i], zh, μ²), 1:num_quark) *
        exp(- PhT² / PhT²avg) /( π * PhT²avg )
end

function get_sf_anselmino2014(data::SidisData)::SidisStructFunc
    f, D = data.f, data.D
    return SidisStructFunc(
        SIDISXSec.zerosf,
        (xB, Q², zh, qT², μ², rtol=0.0) -> FUUT(f, D, xB, zh, qT², μ²),
        SIDISXSec.zerosf,
        SIDISXSec.zerosf,
    )
end

end # module
