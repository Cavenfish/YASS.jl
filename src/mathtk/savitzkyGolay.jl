"""
    savGol(y, ws, order; deriv=0)

Apply a Savitzky-Golay filter to the input signal `y`.

# Arguments
- `y`: Array of data points to be smoothed or differentiated.
- `ws`: Window size (must be odd and ≥ 1).
- `order`: Polynomial order for the filter (must satisfy `ws ≥ order + 2`).
- `deriv`: (Optional) Order of the derivative to compute (default is 0, i.e., smoothing).

# Returns
- Filtered (or differentiated) signal as an array of the same length as `y`.

# Throws
- `ArgumentError` if input arguments are invalid.

# Example
```julia
smoothed = savGol(data, 5, 2)
```
"""
function savGol(y, ws, order; deriv=0)
  #Check inputs
  isodd(ws)      || throw(ArgumentError("window size must be odd"))
  ws ≥ 1         || throw(ArgumentError("window size must be ≥ 1"))
  ws ≥ order + 2 || throw(ArgumentError("ws ≥ order + 2"))
  length(y) > 1  || throw(ArgumentError("len(y) must be > 1"))
  deriv ≥ 0      || throw(ArgumentError("deriv must be ≥ 0"))

  #half window size
  hws = div(ws, 2) 

  #build z vector
  z = collect(-hws:hws)

  #Build Vandermonde matrix
  J = zeros(2*hws+1, order+1)

  #Fill matrix
  for i in 1:order+1
    @. J[:, i] = z ^ (i-1)
  end
  
  #solve for coefficients (adjoint here is the same as transpose)
  #This is the linear least squares solution of J*a = y
  a = inv(J'J) * J' |> (x -> x[deriv+1,:])

  #Pad signal
  b4 = 2*y[1]   .- reverse(y[2:hws+1])
  b5 = 2*y[end] .- reverse(y[end-hws:end-1])
  py = vcat(b4, y, b5)

  #Prep for convolution
  m = length(py)
  n = length(a)
  w = zeros(m + n - 1)

  #Do 1D convolution
  for i in 1:m, j in 1:n
    w[i+j-1] += py[i]*a[j]
  end

  return w[n:end-n+1] #strip padded signal
end