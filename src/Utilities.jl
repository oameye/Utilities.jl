module Utilities

using Reexport
@reexport using DrWatson
using Pipe: @pipe

export @pipe, logspace, geomspace

@doc raw"""
    logspace(start, stop, num=50; base=10, endpoint=true)

Generate a Vector of evenly spaced numbers on a log scale,

```math
v_i = d^{a + (i - 1) \frac{b - a}{n}}
```

for ``i = 1, 2, ⋯, N``, where ``N = `` `num`, ``d = `` `base`, ``a = `` `start`,
``b = `` `stop`, and ``n = `` `num` ``- 1`` if `endpoint` = `true`, otherwise
``n = `` `num`.

# Arguments
- `start::Number`: The starting value of the sequence is `base`^`start`.
- `stop::Number`: The final value of the sequence is `base`^`stop`.
- `num::Integer=50`: Number of samples to generate.

# Keywords
- `base::Number=10`: The base of the log space.
- `endpoint::Bool=true`: If `true`, `stop` is the last sample. Otherwise, it is not included.

# Example
```jldoctest
julia> logspace(-3, 1, 5)
5-element Vector{Float64}:
  0.001
  0.01
  0.1
  1.0
 10.0

```

"""
function logspace(start::Number, stop::Number, num::Integer=50;
                  base::Number=10, endpoint::Bool=true)
    num > 1 || throw(ArgumentError("num <= 1"))
    d = endpoint ? (stop - start) / (num - 1) : (stop - start) / num;
    base = convert(typeof(d), base)
    q = base^d
    res = Vector{typeof(q)}(undef, num)
    res[1] = base^start;
    for i = 2:num-1
        @inbounds res[i] = res[i-1] * q;
    end
    res[num] = endpoint ? base^stop : res[num-1] * q
    return res
end

@doc raw"""
    geomspace(start, stop, num=50; endpoint=true)

Generate a Vector of geometric sequence,

```math
v_i = a \left(\frac{b}{a}\right)^{(i - 1)/n}
```

for ``i = 1, 2, ⋯, N``, where ``N = `` `num`, ``a = `` `start`, ``b = `` `stop`, and
``n = `` `num` ``- 1`` if `endpoint` = `true`, otherwise ``n = `` `num`.

# Arguments
- `start::Number`: The starting value of the sequence.
- `stop::Number`: The final value of the sequence.
- `num::Integer=50`: Number of samples to generate.

# Keywords
- `endpoint::Bool=true`: If `true`, `stop` is the last sample. Otherwise, it is not included.

# Example

```jldoctest
julia> geomspace(1, 1e4, 5)
5-element Vector{Float64}:
     1.0
    10.0
   100.0
  1000.0
 10000.0
```

"""
function geomspace(start::Number, stop::Number, num::Integer=50; endpoint::Bool=true)
    num > 1 || throw(ArgumentError("num <= 1"))
    q = endpoint ? (stop/start)^(1/(num-1)) : (stop/start)^(1/num);
    res_type = typeof(q * oneunit(promote_type(typeof(start), typeof(stop))))
    res = Vector{res_type}(undef, num)
    res[1] = start;
    for i = 2:num-1
        @inbounds res[i] = res[i-1] * q;
    end
    res[num] = endpoint ? stop : res[num-1] * q
    return res
end


end
