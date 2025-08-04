module Constants

export αEM₀
export me, mμ, mτ
export mu, md, ms, mc, mb
export Mp, Md, Mπ, MK, MHe3

const αEM₀ = 1 / 137.036

# masses (in GeV)
# lepton / quark masses with m ; hadron masses with M

const me = 0.000511
const mμ = 0.105658
const mτ = 1.77686

const mu = 0.00216
const md = 0.00467
const ms = 0.0934
const mc = 1.27
const mb = 4.18

const Mp = 0.93827
const Md = 1.87561
const Mπ = 0.13957
const MK = 0.49368
const MHe3 = 3.016 * 0.93149

end # module
