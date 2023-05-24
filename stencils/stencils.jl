module FDStencils

using Symbolics

export GetCoefficient, GetDerivsFD, ReplaceGrids

function GetCoefficient(samples, values)
  @variables c(..)
  npts = length(samples)
  eqns = [sum([c(i)*samples[i]^(j-1) for i in 1:npts]) ~ values[j] for j in 1:npts]
  return Symbolics.solve_for(eqns, [c(i) for i in 1:npts])
end

function GetDerivsFD(sample, order=1)
  @variables gf[1:length(sample)];
  values = [(i == 2) ? 1 : 0 for i in 1:length(sample)]
  coeff = GetCoefficient(sample, values)
  expr = string(sum(coeff[i] * gf[i] for i in 1:length(coeff)))
  for replacerule in ReplaceGrids(length(sample))
      expr = replace(expr, replacerule)
  end
  expr
end

function ReplaceGrids(npts)
    mpt = floor(Int, npts/2) + 1
    # print(mpt)
    [
        if i < mpt
            "gf["*string(i)*"]" => "gf(p.I - "*string(mpt-i)*"*p.DI[dir])"
        else
            "gf["*string(i)*"]" => "gf(p.I + "*string(i-mpt)*"*p.DI[dir])"
        end
    for i in 1:npts]
end

end
