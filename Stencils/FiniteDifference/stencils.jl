module FDStencils

using Symbolics

export GetCoefficient

@variables c(..)

function GetCoefficient(samples, values)
  npts = length(samples)
  eqns = [sum([c(i)*samples[i]^(j-1) for i in 1:npts]) ~ values[j] for j in 1:npts]
  return Symbolics.solve_for(eqns, [c(i) for i in 1:npts])
end

end
