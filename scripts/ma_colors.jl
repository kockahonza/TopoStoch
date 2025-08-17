# figure sizes for publication
single_col_width = 324 * 1.051437737 # corresponds to 511 pdf pt which is ~9cm
double_col_width = single_col_width*2
golden_ratio = 1.618
pt_per_mm_ratio = 3.779527559

# colors to use everywhere
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
