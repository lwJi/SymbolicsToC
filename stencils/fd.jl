module FDStencils

using Symbolics

#export GetCoefficient, GetDerivsFD

#=
GetCoefficient: get stencil for finite difference operator.

  1st order derivs: express f'(0) as the following
    f'(0) = ( c1*f(-2h) + c2*f(-h) + c3*f(h) + c4*f(2h) ) / h,
    which should work for the following base functions: {1, x, x^2, x^3}.
    Hence we have the following linear system to solve for ci:
      ( c1         + c2        + c3       + c4        ) / h = 0,
      ( c1*(-2h)   + c2*(-h)   + c3*(h)   + c4*(2h)   ) / h = (1!),
      ( c1*(-2h)^2 + c2*(-h)^2 + c3*(h)^2 + c4*(2h)^2 ) / h = 0,
      ( c1*(-2h)^3 + c2*(-h)^3 + c3*(h)^3 + c4*(2h)^3 ) / h = 0.

  3rd order derivs: express f'''(0) as the following
    f'''(0) = ( c1*f(-2h) + c2*f(-h) + c3*f(h) + c4*f(2h) ) / h^3,
    which should work for the following base functions: {1, x, x^2, x^3}.
    Hence we have the following linear system to solve for ci:
      ( c1         + c2        + c3       + c4        ) / h^3 = 0,
      ( c1*(-2h)   + c2*(-h)   + c3*(h)   + c4*(2h)   ) / h^3 = 0,
      ( c1*(-2h)^2 + c2*(-h)^2 + c3*(h)^2 + c4*(2h)^2 ) / h^3 = 0,
      ( c1*(-2h)^3 + c2*(-h)^3 + c3*(h)^3 + c4*(2h)^3 ) / h^3 = (3!).
  ...

  Note that h can be cancelled from both sides in the linear system.

=#
function GetDerivsFD(samples, order=1)
  npts = length(samples)
  @variables gf[1:npts], h;
  values = [(i == order+1) ? factorial(order) : 0 for i in 1:npts]
  coeff = GetCoefficient(samples, values)
  return string(simplify(sum(coeff[i] * gf[i] for i in 1:npts))/(h^order))
end

function GetCoefficient(samples, values)
  npts = length(samples)
  #@variables c(..)
  #eqns = [sum([c(i)*samples[i]^(j-1) for i in 1:npts]) ~ values[j]
  #        for j in 1:npts]
  #return Symbolics.solve_for(eqns, [c(i) for i in 1:npts])
  A = transpose(reduce(hcat, [[samples[i]^(j-1) for i in 1:npts]
                              for j in 1:npts]))
  B = transpose(reduce(hcat, values))
  # solve for A*u = B
  return A\B
end

end
