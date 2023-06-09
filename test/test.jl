# test.jl
# (c) Liwe Ji, 05/2023

using Dates

include("../writeio.jl")
using .WriteIO
include("../codes/carpetx.jl")
using .CarpetX

#=
  Wrap functions
=#
function PrintDerivsFunction(pr::Function, funcname::String, sample, order=1)
  pr("template <typename T>")
  pr("inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE T")
  pr(funcname*"(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {")
  pr("  return "*PrintDerivsFD(sample, order)*";")
  pr("}")
  pr()
end

#=
  Write to files
=#
filename = "test.hxx";
thornname = "WaveToyHigherOrderX";

function headpart(pr::Function)
  pr("#ifndef " * replace(uppercase(filename), "."=>"_"))
  pr("#define " * replace(uppercase(filename), "."=>"_"))
  pr()
  pr("/* " * filename * " */")
  pr("/* (c) Liwei Ji " * string(Dates.today()) * " */")
  pr("/* Produced with Julia */")
  pr()
  pr("#include <cctk.h>")
  pr()
  pr("namespace " * thornname * " {")
  pr()
end

function bodypart(pr::Function)
  PrintDerivsFunction(pr, "calc_fd_c1D2O", [-1//1,  0, 1])
  PrintDerivsFunction(pr, "calc_fd_l1D1O", [-1//1,  0])
  PrintDerivsFunction(pr, "calc_fd_r1D1O", [ 0//1,  1])
  PrintDerivsFunction(pr, "calc_fd_c2D2O", [-1//1,  0, 1], 2)
  PrintDerivsFunction(pr, "calc_fd_l2D1O", [-2//1, -1, 0], 2)
  PrintDerivsFunction(pr, "calc_fd_r2D1O", [ 0//1,  1, 2], 2)
  #PrintDerivsFunction(pr, "calc_fd_c1D4O", [-2//1, -1, 0, 1, 2])
  #PrintDerivsFunction(pr, "calc_fd_c2D4O", [-2//1, -1, 0, 1, 2], 2)
end

function tailpart(pr::Function)
  pr("} // namespace " * thornname)
  pr()
  pr("#endif // #ifndef " * replace(uppercase(filename), "."=>"_"))
end

#=
  Main
=#
WriteFile(headpart, bodypart, tailpart, filename)
