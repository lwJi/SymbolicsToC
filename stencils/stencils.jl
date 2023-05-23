module FDStencils

using Symbolics

export GetCoefficient, ReplaceGrids

@variables c(..)

function GetCoefficient(samples, values)
  npts = length(samples)
  eqns = [sum([c(i)*samples[i]^(j-1) for i in 1:npts]) ~ values[j] for j in 1:npts]
  return Symbolics.solve_for(eqns, [c(i) for i in 1:npts])
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
