#= Part of SIDISStructFunc/Coll.jl                                                                =#

#= Example usage:
using SIDISXSec, SIDISXSec.Constants
data = get_sidis_data(Mp, "NNPDF31_nnlo_as_0118_1000", Mπ, "DSS14_pion_plus_nlo", A=1, Z=1)
FUUT = SIDISXSec._get_sf_coll(data).FUUT
SIDISXSec.Coll.generate_coll(FUUT, "FUUT_lo", "proton", "pi+", "NNPDF31_nnlo_as_0118_1000,DSS14_pion_plus_nlo")
=#

const __rtol = 1e-4

# manually defined bounds
const xB_low, zh_low = 0.001, 0.1
const qT²oQ²_low, qT²oQ²_up′ = 0.01^2, 10.0^2
const μ²_low, μ²_up = 1.0^2, 100.0^2
# deduced bounds
const qT²oQ²_up = get_qT²oQ²max(xB_low, zh_low)

grid_decode(qT²oQ², R, A) = let
    qT²oQ²_max = qT²oQ²^(R) * qT²oQ²_up^(1-R)
    xB = inv( 1 + zh_low/(1-zh_low) * qT²oQ²_max^(1-A) * qT²oQ²_up^(  A) )
    zh = inv( 1 + xB_low/(1-xB_low) * qT²oQ²_max^(  A) * qT²oQ²_up^(1-A) )
    return xB, zh, qT²oQ²
end
grid_encode(xB, zh, qT²oQ²) = let
    qT²oQ²_max = get_qT²oQ²max(xB,zh)
    R = log( qT²oQ²_up / qT²oQ²     , qT²oQ²_up / qT²oQ²_max               )
    A = log( qT²oQ²_up / qT²oQ²_max , qT²oQ²_up / get_qT²oQ²max(xB_low,zh) )
    return qT²oQ², R, A
end

shift_scale(μ², Cμ, μ_cut=1.0) = max(μ_cut, Cμ^2 * μ²)

function generate_coll(sf::Function, sfname::String, Nname::String, hname::String, paramname::String; Cμ=1.0, rtol=__rtol)
    setname = "Coll_$(sfname)_$(Nname)_$(hname)"
    desc = "Coll $sfname with initial state $Nname and final state $hname, using parametrization $paramname. [rtol=$(@sprintf("%.1E", rtol))] [μ=$Cμ*Q] [ qT²oQ²_max = qT²oQ²^(R) * qT²oQ²_up^(1-R); xB = inv( 1 + zh_low/(1-zh_low) * qT²oQ²_max^(1-A) * qT²oQ²_up^(  A) ); zh = inv( 1 + xB_low/(1-xB_low) * qT²oQ²_max^(  A) * qT²oQ²_up^(1-A) ); qT²oQ²_up=get_qT²oQ²max(xB_low, zh_low); xB_low=$xB_low; zh_low=$zh_low ]"

    print("\n")
    println("Generating $setname")
    flush(stdout)

    # Test -----------------------------------------------------------------------------------------

    println("Test: $sfname(0.3, 10, 0.3, 1, 10, 5e-2) = $(sf(0.3, 10, 0.3, 1, 10, 5e-2))")
    flush(stdout)

    # Calculations ---------------------------------------------------------------------------------
    qT²oQ²_grid = logrange(qT²oQ²_low, qT²oQ²_up′, length=50+1)
    R_grid = range(0, 1, length=50+1)
    A_grid = range(0, 1, length=80+1)
    μ²_grid = logrange(μ²_low, μ²_up, length=30+1)
    values = Array{TMDGridFloat}(undef, length(qT²oQ²_grid), length(R_grid), length(A_grid), length(μ²_grid))

    println("Number of threads: $(Threads.nthreads())")
    print("\r0% done")
    flush(stdout)
    start_time = time_ns()
    μ²_count = Threads.Atomic{Int}(0)

    Threads.@threads for iμ² ∈ eachindex(μ²_grid)
        for iA ∈ eachindex(A_grid), iqT²oQ² ∈ eachindex(qT²oQ²_grid)
            # R = 1 (qT²oQ² = qT²oQ²_max)
            values[iqT²oQ², end, iA, iμ²] = 0
            # R = 0 (qT²oQ²_max = qT²oQ²_up)
            μ² = μ²_grid[iμ²]
            value = sf(xB_low, μ², zh_low, qT²oQ²_grid[iqT²oQ²]*μ², shift_scale(μ², Cμ), rtol)
            values[iqT²oQ², begin, iA, iμ²] = value
            if isnan(value) @warn "Encountered NaN." end
        end
        for iA ∈ eachindex(A_grid), iR ∈ 2:lastindex(R_grid)-1, iqT²oQ² ∈ eachindex(qT²oQ²_grid)
            xB, zh, qT²oQ² = grid_decode(qT²oQ²_grid[iqT²oQ²], R_grid[iR], A_grid[iA])
            μ² = μ²_grid[iμ²]
            value = sf(xB, μ², zh, qT²oQ²*μ², shift_scale(μ², Cμ), rtol)
            values[iqT²oQ², iR, iA, iμ²] = value
            if isnan(value) @warn "Encountered NaN." end
        end
        Threads.atomic_add!(μ²_count, 1)
        print("\r$(floor(Int, 100μ²_count[]/lastindex(μ²_grid)))% done")
        flush(stdout)
    end

    end_time = time_ns()
    println(". $(ceil(Int, (end_time-start_time)/1e9))s used")
    flush(stdout)

    # Write files ----------------------------------------------------------------------------------
    write_tmdgrid(setname,
        TMDGrid([
            (name="qT²oQ²", type="log", start=qT²oQ²_low, stop=qT²oQ²_up′, length=length(qT²oQ²_grid)),
            (name="R",      type="lin", start=0.0,        stop=1.0,        length=length(R_grid)     ),
            (name="A",      type="lin", start=0.0,        stop=1.0,        length=length(A_grid)     ),
            (name="μ²",     type="log", start=μ²_low,     stop=μ²_up,      length=length(μ²_grid)    ),
            ], values, desc)
    )
end
function generate_coll(no::Int, fname::String, Nname::String, Dname::String, hname::String)
    generate_coll(no, fname, Nname, 1, 1, Dname, hname)
end
