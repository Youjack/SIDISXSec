module QCDData

using SpecialFunctions: beta
using LaTeXStrings
using cLHAPDF

export num_quark, quark_name, quark_charge, quark_code, quark_code2idx
export iso_code, iso_list
export get_pdf_toy
export get_qcd_data, get_αs, get_qcd_density

const num_quark = 10
const quark_name     =     ( L"u"      , L"\bar{u}", L"d"      , L"\bar{d}", L"s"      , L"\bar{s}", L"c"      , L"\bar{c}" , L"b"      , L"\bar{b}" )
const quark_charge   =     ( +2/3      , -2/3      , -1/3      , +1/3      , -1/3      , +1/3      , +2/3      , -2/3       , -1/3      , +1/3       )
const quark_code     =     (  2        , -2        ,  1        , -1        ,  3        , -3        ,  4        , -4         ,  5        , -5         )
const quark_code2idx = Dict(  2 => 1   , -2 => 2   ,  1 => 3   , -1 => 4   ,  3 => 5   , -3 => 6   ,  4 => 7   , -4 => 8    ,  5 => 9   , -5 => 10   )

"Isospin transform of mccode."
iso_code(id) =
    if     id ==  2 return  1
    elseif id ==  1 return  2
    elseif id == -2 return -1
    elseif id == -1 return -2
    else return id end
"Isospin transform of list of quarks."
iso_list(val) = [ val[i] for i =
    [ 3, 4, 1, 2, 5, 6, 7, 8, 9, 10 ]
]

"""
Get toy PDF. \\
`x` is not constrained to be in `(0,1)`.
"""
function get_pdf_toy(id::Integer)::Function
    av = 0.5; bv = 3.0
    as = -0.08; bs = 7.0
    xv(x) = 1/x * x^av * (1-x)^bv / beta(av, bv+1)
    xs(x) = 1/x * x^as * (1-x)^bs / beta(as+1, bs+1)
    if id == 2 # u
        return x -> 2.0xv(x) + 0.03xs(x)
    elseif id == 1 # d
        return x -> 1.0xv(x) + 0.036xs(x)
    elseif id == -2 # ū
        return x -> 0.03xs(x)
    elseif id == -1 # d̄
        return x -> 0.036xs(x)
    elseif id == 4 || id == -4 # c c̄
        return x -> 0.005xs(x)
    elseif id == 3 || id == -3 # s s̄
        return x -> 0.016xs(x)
    else
        return x -> 0.0
    end
end

#= qcd_data from LHAPDF ===========================================================================#

"""
    get_qcd_data(setname_nmem::String)::Ptr{Cvoid}
"""
function get_qcd_data(setname_nmem::String)::Ptr{Cvoid}
    mkPDF(setname_nmem)
end

"""
    get_αs(qcd_data::Ptr{Cvoid}, μ²)::Float64

Get αs at scale `μ²`.
"""
function get_αs(qcd_data::Ptr{Cvoid}, μ²::Real)::Float64
    alphasQ2(qcd_data, μ²)
end

"""
    get_qcd_density(qcd_data::Ptr{Cvoid}, id, μ²,
        conj=false, A=1, Z=1)::Function

Get PDF or FF at scale `μ²`, as f(x).
(`x` is not constrained to be in `(0,1)`)
"""
function get_qcd_density(qcd_data::Ptr{Cvoid}, id::Integer, μ²::Real,
        conj=false, A=1, Z=1)::Function
    x -> 1/x * (
        +    Z  * xfxQ2(qcd_data, (-1)^conj *          id,  x, μ²)
        + (A-Z) * xfxQ2(qcd_data, (-1)^conj * iso_code(id), x, μ²)
        )/A
end

end # module
