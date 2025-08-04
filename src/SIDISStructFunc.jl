#= Part of SIDISXsec.jl                                                                           =#
#= SIDIS structure functions                                                                      =#
# Following conventions in [bacchetta2007semiinclusive]

"""
SIDIS structure functions. \\
Unneeded structure functions can be set to be `zeroSF`.
- `F(xB, Q², zh, qT², μ², rtol=_rtol)`
"""
struct SidisStructFunc{F1,F2,F3,F4}
    FUUL      :: F1
    FUUT      :: F2
    FUUcosϕh  :: F3
    FUUcos2ϕh :: F4
end
Base.show(io::IO, SF::SidisStructFunc) = print(io, "SidisStructFunc{...}")

function zerosf(xB, Q², zh, qT², μ², rtol=_rtol)::Float64
    0.0
end
get_sf_zero() = SidisStructFunc(zerosf, zerosf, zerosf, zerosf)

#= TMD factorization ==============================================================================#

include("SIDISStructFunc/Anselmino2014.jl")
using .Anselmino2014

include("SIDISStructFunc/Barone2015.jl")
using .Barone2015

include("SIDISStructFunc/Bastami2019.jl")
using .Bastami2019

include("SIDISStructFunc/SF_SV19.jl")
using .SF_SV19

# include("SIDISStructFunc/SF_PV17.jl")
# using .SF_PV17

include("SIDISStructFunc/TMD.jl")
using .TMD

#= Collinear factorization ========================================================================#

include("SIDISStructFunc/Coll.jl")
using .Coll

#= Matching prescription ==========================================================================#

include("SIDISStructFunc/WY.jl")
using .WY

include("SIDISStructFunc/InEW.jl")
using .InEW
