#= Part of SIDISStructFunc/TMD.jl                                                                 =#

#= Example usage:
using SIDISXSec, SIDISXSec.Constants
data = get_sidis_data(Mp, "NNPDF31_nnlo_as_0118_1000_splitext", Mπ, "DSS14_pion_plus_nlo_splitext", A=1, Z=1)
FUUT = SIDISXSec.get_sf_sv19(data).FUUT
SIDISXSec.TMD.generate_tmd(FUUT, "FUUT_lo", "proton", "pi+", "SV19:NNPDF31_nnlo_as_0118_1000,DSS14_pion_plus_nlo")
=#

const xB_low, zh_low = 0.001, 0.1
const qT²_low = (1e-3)^2
const Q²_low, Q²_up = 1.0^2, 100.0^2

const __rtol = 1e-4

const qT²oQ²_up = 2^2
grid_decode(R, Q²) = let
    qT² = qT²_low^(1-R) * (qT²oQ²_up*Q²)^R
    return Q², qT²
end
grid_encode(Q², qT²) = let
    R = log( qT²oQ²_up*Q² / qT²_low, qT² / qT²_low )
    return R, Q²
end

function generate_tmd(sf::Function, sfname::String, Nname::String, hname::String, paramname::String;
        path=nothing, rtol=__rtol)
    setname = "TMD_$(sfname)_$(Nname)_$(hname)"
    desc = "TMD $sfname with initial state $Nname and final state $hname, using parametrization $paramname. [rtol=$(@sprintf("%.1E", rtol))] [μ²=Q²] [ qT² = qT²_low^(1-R) * (qT²oQ²_up*Q²)^R; qT²_low=$qT²_low; qT²oQ²_up=$qT²oQ²_up ]"

    print("\n")
    println("Generating $setname")
    flush(stdout)

    # Test -----------------------------------------------------------------------------------------

    println("Test: $sfname(0.3, 10, 0.3, 1, 10, 5e-2) = $(sf(0.3, 10, 0.3, 1, 10, 5e-2))")
    flush(stdout)

    # Calculations ---------------------------------------------------------------------------------
    R_grid = range(0.0, 1.0, length=60+1)
    xB_grid = logrange(xB_low, 1.0, length=80+1)
    zh_grid = logrange(zh_low, 1.0, length=40+1)
    Q²_grid = logrange(Q²_low, Q²_up, length=30+1)
    values = Array{TMDGridFloat}(undef, length(R_grid), length(xB_grid), length(zh_grid), length(Q²_grid))

    println("Number of threads: $(Threads.nthreads())")
    print("\r0% done")
    flush(stdout)
    start_time = time_ns()
    Q²_count = Threads.Atomic{Int}(0)

    Threads.@threads for iQ² ∈ eachindex(Q²_grid)
        for izh ∈ eachindex(zh_grid), ixB ∈ eachindex(xB_grid), iR ∈ eachindex(R_grid)
            Q², qT² = grid_decode(R_grid[iR], Q²_grid[iQ²])
            values[iR, ixB, izh, iQ²] = sf(xB_grid[ixB], Q², zh_grid[izh], qT², Q², rtol)
        end
        Threads.atomic_add!(Q²_count, 1)
        print("\r$(floor(Int, 100Q²_count[]/lastindex(Q²_grid)))% done")
        flush(stdout)
    end

    end_time = time_ns()
    println(". $(ceil(Int, (end_time-start_time)/1e9))s used")
    flush(stdout)

    # Write files ----------------------------------------------------------------------------------
    write_tmdgrid(setname,
        TMDGrid([
            (name="R",  type="lin", start=0.0,    stop=1.0,   length=length(R_grid) ),
            (name="xB", type="log", start=xB_low, stop=1.0,   length=length(xB_grid)),
            (name="zh", type="log", start=zh_low, stop=1.0,   length=length(zh_grid)),
            (name="Q²", type="log", start=Q²_low, stop=Q²_up, length=length(Q²_grid)),
            ], values, desc),
        path=path
    )
end
