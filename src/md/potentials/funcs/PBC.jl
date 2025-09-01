
function pbcVec(u::AbstractVector, i::Int64, j::Int64, lat::AbstractMatrix)
  r  = u[j] .- u[i]
  d  = norm(r)

  if iszero(lat)
    return r,d
  end

  v  = zeros(3)
  rp = zeros(3)
  rn = zeros(3)

  for k = 1:3
    v  .= lat[k,:]
    rp .= (u[j] .+ v) .- u[i]
    rn .= (u[j] .- v) .- u[i]
    dp  = norm(rp)
    dn  = norm(rn)

    if dp < d && dp < dn
      d  = dp
      r .= rp
    elseif dn < d
      d  = dn
      r .= rn
    end
  end

  r, d
end

function _pbc!(
  F::AbstractVector, u::AbstractVector, a::Int64, b::Int64,
  func::FT, L::AbstractMatrix, NC::Vector{Int64}, p::PT
) where {FT, PT}
  E = 0.0
  for i = -NC[1]:NC[1]
    for j = -NC[2]:NC[2]
      for k = -NC[3]:NC[3]
        (i,j,k) == (0,0,0) && continue

        r2     = u[b] + (L[1, :] * i) + (L[2, :] * j) + (L[3, :] * k)

        e,f    = func(u[a], r2, p...)
        E     += e
        F[a] .-= f
        F[b] .+= f
      end
    end
  end
  E
end
