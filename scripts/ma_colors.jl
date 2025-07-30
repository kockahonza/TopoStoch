module MAColors

using Colors
using ColorSchemes

gray1() = colorant"snow3"
lightgray_kabirs() = colorant"#c8c8c8"

ligand_1s() = colorant"#d95f02"

ac_ss() = colorant"#0070b0"
ac_nss_eq() = colorant"#e69e00"
ac_nss_neq() = colorant"#009c73"

bs_cmap() = ColorSchemes.RdBu_4
bs_0s() = ColorSchemes.RdBu_4[end]
bs_1s() = ColorSchemes.RdBu_4[begin]
bs_else() = colorant"navajowhite"

end
