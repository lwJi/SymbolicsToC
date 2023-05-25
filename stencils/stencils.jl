module FDStencils

using Symbolics

export GetCoefficient, GetDerivsFD, ReplaceGrids

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
function GetDerivsFD(sample, order=1)
  @variables gf[1:length(sample)], h;
  values = [(i == order+1) ? factorial(order) : 0 for i in 1:length(sample)]
  coeff = GetCoefficient(sample, values)
  expr = string(simplify(sum(coeff[i] * gf[i] for i in 1:length(coeff)))/(h^order))
  for replacerule in ReplaceGrids(length(sample))
      expr = replace(expr, replacerule)
  end
  expr
end

function GetCoefficient(samples, values)
  @variables c(..)
  npts = length(samples)
  eqns = [sum([c(i)*samples[i]^(j-1) for i in 1:npts]) ~ values[j] for j in 1:npts]
  return Symbolics.solve_for(eqns, [c(i) for i in 1:npts])
end

#=
For CarpetX
=#
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
