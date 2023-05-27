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
function GetDerivsFD(samples, order=1)
  npts = length(samples)
  @variables gf[1:npts], h;
  values = [(i == order+1) ? factorial(order) : 0 for i in 1:npts]
  coeff = GetCoefficient(samples, values)
  expr = string(simplify(sum(coeff[i] * gf[i] for i in 1:npts))/(h^order))
  pt0 = findall(x->x == 0, samples)[1]
  for replacerule in ReplaceGrids(npts, pt0)
    expr = replace(expr, replacerule)
  end
  expr
end

function GetCoefficient(samples, values)
  npts = length(samples)
  #@variables c(..)
  #eqns = [sum([c(i)*samples[i]^(j-1) for i in 1:npts]) ~ values[j] for j in 1:npts]
  #return Symbolics.solve_for(eqns, [c(i) for i in 1:npts])
  A = transpose(reduce(hcat, [[samples[i]^(j-1) for i in 1:npts] for j in 1:npts]))
  B = transpose(reduce(hcat, values))
  # solve for A*u = B
  A\B
end

#=
For CarpetX
=#
function ReplaceGrids(npts, pt0)
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

end
