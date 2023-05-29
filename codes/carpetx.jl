module CarpetX

include("../stencils/fd.jl")
using .FDStencils

export PrintDerivsFD

function ReplaceGFNames(npts, pt0)
  # print(pt0)
  [ if i == pt0
      "gf["*string(i)*"]" => "gf(p.I)"
    elseif abs(i-pt0) == 1
      if i < pt0
        "gf["*string(i)*"]" => "gf(p.I - p.DI[dir])"
      else
        "gf["*string(i)*"]" => "gf(p.I + p.DI[dir])"
      end
    elseif i < pt0
      "gf["*string(i)*"]" => "gf(p.I - "*string(pt0-i)*"*p.DI[dir])"
    else
      "gf["*string(i)*"]" => "gf(p.I + "*string(i-pt0)*"*p.DI[dir])"
    end
  for i in 1:npts]
end

function PrintDerivsFD(samples, order=1)
  expr = FDStencils.GetDerivsFD(samples, order)
  pt0 = findall(x->x == 0, samples)[1]
  for replacerule in ReplaceGFNames(length(samples), pt0)
    expr = replace(expr, replacerule)
  end
  return expr
end

end
