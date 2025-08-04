module QCDData

using SpecialFunctions: beta
using LaTeXStrings, Printf
using cLHAPDF
using Interpolations

export num_quark, quark_name, quark_charge, quark_code, quark_code2idx
export iso_code, iso_list
export get_pdf_toy
export get_qcd_data, get_αs, get_qcd_density

export TMDGrid, write_tmdgrid, read_tmdgrid, interpolate_tmdgrid

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
    get_qcd_density(qcd_data::Ptr{Cvoid}, id, μ²)::Function

Get PDF or FF at scale `μ²`, as f(x). \\
`x` is not constrained to be in `(0,1)`.
"""
function get_qcd_density(qcd_data::Ptr{Cvoid}, id::Integer, μ²::Real)::Function
    x -> 1/x * xfxQ2(qcd_data, id, x, μ²)
    #= x -> 0<x && x<1 ?
        1/x * xfxQ2(qcd_data, id, x, μ²) :
        throw(DomainError(x, "momentum fraction out of range (QCD)")) =#
end

#= TMDGrid ========================================================================================#

const TMDGRID_PATH = split(ENV["TMDGRID_PATH"], ":")[1]
const VarGrid = @NamedTuple{name::String, type::String, start::Float64, stop::Float64, length::Int}

"""
Grid to store TMD data.
- `vars = [(name::String, type::String, start::Float64, stop::Float64, length::Int), ...]`
- `values` is a multi-dimensional array. It mimics the function `values(var1, var2, ...)`.
"""
struct TMDGrid
    vars :: Vector{VarGrid}
    values
    function TMDGrid(vars, values::AbstractArray)
        size_values = size(values)
        if length(vars) != length(size_values)
            throw(ErrorException("length(vars) != length(size(values))"))
        end
        for n ∈ eachindex(vars)
            if vars[n].length != size_values[n]
                throw(ErrorException("vars[n].length != size(values)[n]"))
            end
        end
        new(vars, values)
    end
end
Base.show(io::IO, tmdgrid::TMDGrid) = print(io,
    "TMDGrid with vars: ", join([var_grid.name for var_grid ∈ tmdgrid.vars], ", ")
)

const _sigdigits = 8
round2str(x, sigdigits) = Printf.format(Printf.Format("%$(sigdigits+8).$(sigdigits-1)E"), x)

function write_tmdgrid(name::String, tmdgrid::TMDGrid; desc="No description.", sigdigits=_sigdigits)
    open(joinpath(TMDGRID_PATH, name*".tmdgrid"), "w") do io
        # write description
        write(io, desc*"\n")
        # write vars
        for var_grid ∈ tmdgrid.vars
            write(io, rpad(var_grid.name, 10))
            write(io, ":")
            write(io, lpad(var_grid.type, 10))
            write(io, round2str(var_grid.start, sigdigits))
            write(io, round2str(var_grid.stop,  sigdigits))
            write(io, lpad("$(var_grid.length)", 5))
            write(io, "\n")
        end
        # write values
        col = 1
        for I ∈ CartesianIndices(tmdgrid.values)
            if col != Tuple(I)[end]
                write(io, "\n")
                col += 1
            end
            write(io, round2str(tmdgrid.values[I], sigdigits))
        end
        write(io, "\n")
    end
    return nothing
end

function read_tmdgrid(name::String)::TMDGrid
    open(joinpath(TMDGRID_PATH, name*".tmdgrid"), "r") do io
        readline(io) # skip description
        split_line::Vector{SubString{String}} = Vector{SubString{String}}(undef, 0)
        # read vars
        vars = Vector{VarGrid}(undef, 0)
        var_lengths = Vector{Int}(undef, 0)
        while (split_line = split(readline(io)))[2] == ":"
            var_grid = (
                name = split_line[1], type = split_line[3],
                start = parse(Float64, split_line[4]), stop = parse(Float64, split_line[5]),
                length = parse(Int, split_line[6]))
            push!(vars, var_grid)
            push!(var_lengths, var_grid.length)
        end
        # read values
        values_size = Tuple(var_lengths)
        values_ndims = length(values_size)
        values = Array{Float64,values_ndims}(undef, values_size)
        line_size = Tuple(@view var_lengths[1:end-1])
        for col ∈ 1:values_size[end]
            selectdim(values, values_ndims, col) .=
                reshape(parse.(Float64, split_line), line_size)
            split_line = split(readline(io))
        end
        return TMDGrid(vars, values)
    end
end

function get_node(x, var_grid::VarGrid)::Float64
    xl, xu, L = var_grid.start, var_grid.stop, var_grid.length
    if     var_grid.type == "lin"
        return round(( (xu-x) + L * (x-xl) )/ (xu-xl), sigdigits=8)
    elseif var_grid.type == "log"
        return round(( log(xu/x) + L * log(x/xl) )/ log(xu/xl), sigdigits=8)
    else
        throw(ErrorException("get_node: unkown VarGrid type `$(var_grid.type)`."))
    end
end
function interpolate_tmdgrid(tmdgrid::TMDGrid)::Function
    itp = interpolate(tmdgrid.values, BSpline(Cubic(Line(OnGrid()))))
    return (xs...) -> itp((get_node(xs[n], tmdgrid.vars[n]) for n ∈ eachindex(xs))...)
end

end # module
