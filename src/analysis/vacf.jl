"""
    vacfInps

Input structure for velocity autocorrelation function (VACF) calculations.

# Fields
- `vel`: Velocity data.
- `mas`: Masses.
- `Hz`: Sampling frequency.
- `norm`: Normalize flag.
- `win`: Window function.
- `pad`: Padding factor.
- `mir`: Mirror flag.
"""
struct vacfInps{V,B,W, I<:Integer, F<:AbstractFloat}
  vel::V
  mas::Vector{F}
  Hz::F
  norm::B
  win::W
  pad::I
  mir::B
end

"""
    vacfOut

Output structure for VACF and VDOS calculations.

# Fields
- `c`: Raw VACF.
- `C`: Windowed VACF.
- `v`: Frequency axis.
- `I`: Intensity.
"""
struct vacfOut{F<:AbstractFloat}
  c::Vector{F}
  C::Vector{F}
  v::Vector{F}
  I::Vector{F}
end

Hann( i,N) = sin((pi*i)/N)^2
Welch(i,N) = 1 - ((i-N/2)/(N/2))^2
HannM(i,N) = cos((pi*i)/(2*(N-1)))^2

"""
    window!(out, W)

Apply a window function to the VACF output.

# Arguments
- `out`: `vacfOut` object.
- `W`: Window function.
"""
function window!(out, W)
  N = length(out.C)
  w = [W(i,N) for i in 0:N-1]

  out.C .*= w
end

"""
    mirror!(out)

Mirror the VACF output for improved frequency resolution.

# Arguments
- `out`: `vacfOut` object.
"""
function mirror!(out)
  N::Int = length(out.C) / 2 |> ceil

  out.C .= [reverse(out.C[1:N]); out.C[2:N]]
end

"""
    vacf!(inp, out; atms=nothing)

Compute the velocity autocorrelation function (VACF).

# Arguments
- `inp`: `vacfInps` object.
- `out`: `vacfOut` object.
- `atms`: (Optional) Indices of atoms to include.
"""
function vacf!(inp, out; atms=nothing)
  T = length(inp.vel)       # Total timesteps 
  N = length(inp.vel[1])    # Total number of atoms
  D = length(inp.vel[1][1]) # Total number of dimensions

  if atms == nothing
    atms = 1:N
  end
  
  for i in atms
    for j in 1:D
      data = [a[i][j] for a in inp.vel]
      forw = fft(data)
      tmp  = forw .* conj(forw)
      back = ifft(tmp) ./ T
      
      @views out.c  .+= inp.mas[i] * real(back[1:T-1])
    end
  end

  if inp.norm
    out.C[1:T-1] .= out.c ./ out.c[1]
  else
    out.C[1:T-1] .= out.c
  end

end

"""
    VDOS(inp; atms=nothing)

Compute the vibrational density of states (VDOS) from VACF input.

# Arguments
- `inp`: `vacfInps` object.
- `atms`: (Optional) Indices of atoms to include.

# Returns
- `vacfOut` object with VDOS data.
"""
function VDOS(inp; atms=nothing)
  # Initialize output
  k   = 29979245800.0 # Divide by this to convert from Hz to cm^-1
  c   = length(inp.vel) - 1 |> zeros
  C   = if inp.mir
    (length(c)*inp.pad*2) - 1 |> zeros
  else
    length(c)*inp.pad |> zeros
  end
  N   = length(C)
  tmp = fftfreq(N, inp.Hz) ./ k
  n   = div(length(tmp), 2)
  v   = abs.(tmp[1:n])
  I   = length(v) |> zeros
  out = vacfOut(c,C,v,I)

  vacf!(inp, out; atms=atms)
  window!(out, inp.win)
  
  if inp.mir
    mirror!(out)
  end

  tmp    = fft(out.C)
  out.I .= abs.(tmp[1:n])

  out
end

"""
    getDiffusionCoefficient(out; D=3)

Calculate the diffusion coefficient from VACF output.

# Arguments
- `out`: `vacfOut` object.
- `D`: Dimensionality (default: 3).

# Returns
- Diffusion coefficient (Float64).
"""
function getDiffusionCoefficient(out; D=3)
  N = length(out.C)

  #By default the diffusion coefficient is in units of
  # Angstrom * sqrt(eV / amu) | 0.0098226 is conversion to cm^2/s
  ( sum(out.c) / (D*N) ) * 0.00982269475
end