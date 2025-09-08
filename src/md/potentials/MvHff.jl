"""
CO-CO Potential

Based on:
van Hemert, Marc C., Junko Takahashi, and Ewine F. van Dishoeck. "Molecular
dynamics study of the photodesorption of CO ice." The Journal of Physical
Chemistry A 119.24 (2015): 6354-6369.

Link:
https://pubs.acs.org/doi/full/10.1021/acs.jpca.5b02611
"""
MvHff(; constraints=nothing) = Calculator(MvHff; EF=MvHff!, constraints=constraints)

struct _MvHff_PotVars{F<:Float64} <: PotVars
  D::F
  a::F
  req::F
  Acc::F
  Aco::F
  Aoo::F
  Bcc::F
  Bco::F
  Boo::F
  Ccc::F
  Cco::F
  Coo::F
  Qc::F
  Qo::F
  αc::F
  αo::F
end

MvHff(bdys::Vector{MyAtoms}) = _MvHff_PotVars(
  11.230139012256362,
  2.3281,
  1.1282058129221093,
  361.367206403597,
  1516.76265699823,
  6370.185468304371,
  2.8345891887553925,
  3.5432364859442407,
  4.2518837831330885,
  33.44955570065988,
  15.189133724736946,
  10.546349734130885,
  -1.7835026375774934,
  -2.333732174702465,
  3.843702939952312,
  2.131611069944055
)

function MvHff!(F, u, p)
  E = 0.0
  P = p.potVars

  for mol in p.mols
    E += _Morse!(F, u, mol[1], mol[2], P.D, P.a, P.req)
  end

  for i = 1:length(p.mols)
    c1, o1 = p.mols[i]

    for j = i+1:length(p.mols)
      c2, o2 = p.mols[j]

      #C--C
      E += _Buckingham!(F, u, c1, c2, P.Acc, P.Bcc, P.Ccc)

      #C--O
      E += _Buckingham!(F, u, c1, o2, P.Aco, P.Bco, P.Cco)

      #O--C
      E += _Buckingham!(F, u, o1, c2, P.Aco, P.Bco, P.Cco)

      #O--O
      E += _Buckingham!(F, u, o1, o2, P.Aoo, P.Boo, P.Coo)

      #Special Electrostatics
      E += _electroMvH!(F, u, p.mols[[i,j]], P.Qc, P.Qo, P.αc, P.αo, P.req)
    end
  end

  E
end

function _electroMvH!(F, u, par, Qc, Qo, αc, αo, req)

  #Atom indicies
  C1, O1 = par[1]
  C2, O2 = par[2]

  #Atom locations
  c1, o1 = u[par[1]]
  c2, o2 = u[par[2]]

  #Bond lengths
  R1 = o1 - c1
  R2 = o2 - c2
  r1 = norm(R1)
  r2 = norm(R2)

  #CM weights
  wo = 0.57135
  wc = 0.42865

  #X locations
  x1 = @. (wc*c1) + (wo*o1)
  x2 = @. (wc*c2) + (wo*o2)

  #Variable charges
  Qc1 = Qc * exp(-αc*(r1 - req))
  Qo1 = Qo * exp(-αo*(r1 - req))
  Qx1 = -(Qc1 + Qo1)
  Qc2 = Qc * exp(-αc*(r2 - req))
  Qo2 = Qo * exp(-αo*(r2 - req))
  Qx2 = -(Qc2 + Qo2)

  #The additional contribution to the forces resulting from
  #the position dependent charges is: -α * Q * rvec / r

  #C--C
  ϵ, f   = _Coulomb(c1, c2, Qc1, Qc2)
  E      = ϵ
  F_Q1   = @. αc * ϵ * R1 / r1
  F_Q2   = @. αc * ϵ * R2 / r2
  @. F[C1] -= (f + F_Q1)
  @. F[O1] += F_Q1
  @. F[C2] +=  f - F_Q2
  @. F[O2] += F_Q2

# C-X
  ϵ, f   = _Coulomb(c1, x2, Qc1, Qx2)
  E     += ϵ
  F_Q1   = @. αc * ϵ * R1 / r1
  F_Q2   = @. - (αc * Qc2 + αo * Qo2) * ϵ/Qx2 * R2 / r2
  @. F[C1] -= (f + F_Q1)
  @. F[O1] += F_Q1
  @. F[C2] += (wc * f) - F_Q2
  @. F[O2] += (wo * f) + F_Q2

  # C-O
  ϵ, f   = _Coulomb(c1, o2, Qc1, Qo2)
  E     += ϵ
  F_Q1   = @. αc * ϵ * R1 / r1
  F_Q2   = @. αo * ϵ * R2 / r2
  @. F[C1] -= (f + F_Q1)
  @. F[O1] += F_Q1
  @. F[C2] -= F_Q2
  @. F[O2] += (f + F_Q2)

  # X-C
  ϵ, f   = _Coulomb(x1, c2, Qx1, Qc2)
  E     += ϵ
  F_Q1   = @. - (αc * Qc1 + αo * Qo1) * ϵ/Qx1 * R1 / r1
  F_Q2   = @. αc * ϵ * R2 / r2
  @. F[C1] -= ((wc * f) + F_Q1)
  @. F[O1] += F_Q1 - (wo * f)
  @. F[C2] += (f - F_Q2)
  @. F[O2] += F_Q2

  # X-X
  ϵ, f   = _Coulomb(x1, x2, Qx1, Qx2)
  E     += ϵ
  F_Q1   = @. - (αc * Qc1 + αo * Qo1) * ϵ/Qx1 * R1 / r1
  F_Q2   = @. - (αc * Qc2 + αo * Qo2) * ϵ/Qx2 * R2 / r2
  @. F[C1] += ((-wc * f) - F_Q1)
  @. F[O1] += ((-wo * f) + F_Q1)
  @. F[C2] += (( wc * f) - F_Q2)
  @. F[O2] += (( wo * f) + F_Q2)


  # X-O
  ϵ, f   = _Coulomb(x1, o2, Qx1, Qo2)
  E     += ϵ
  F_Q1   = @. - (αc * Qc1 + αo * Qo1) * ϵ/Qx1 * R1 / r1
  F_Q2   = @. αo * ϵ * R2 / r2
  @. F[C1] -= ((wc * f) + F_Q1)
  @. F[O1] += (F_Q1 - (wo * f))
  @. F[C2] -= F_Q2
  @. F[O2] += (f + F_Q2)

  # O-C
  ϵ, f   = _Coulomb(o1, c2, Qo1, Qc2)
  E     += ϵ
  F_Q1   = @. αo * ϵ * R1 / r1
  F_Q2   = @. αc * ϵ * R2 / r2
  @. F[C1] -= F_Q1
  @. F[O1] += (F_Q1 - f)
  @. F[C2] += (f - F_Q2)
  @. F[O2] += F_Q2

  # O-X
  ϵ, f   = _Coulomb(o1, x2, Qo1, Qx2)
  E     += ϵ
  F_Q1   = @. αo * ϵ * R1 / r1
  F_Q2   = @. - (αc * Qc2 + αo * Qo2) * ϵ/Qx2 * R2 / r2
  @. F[C1] -= F_Q1
  @. F[O1] += (F_Q1 - f)
  @. F[C2] += ((wc * f) - F_Q2)
  @. F[O2] += ((wo * f) + F_Q2)

  # O-O
  ϵ, f   = _Coulomb(o1, o2, Qo1, Qo2)
  E     += ϵ
  F_Q1   = @. αo * ϵ * R1 / r1
  F_Q2   = @. αo * ϵ * R2 / r2
  @. F[C1] -= F_Q1
  @. F[O1] += (F_Q1 - f)
  @. F[C2] -= F_Q2
  @. F[O2] += (f + F_Q2)
  

  E
end