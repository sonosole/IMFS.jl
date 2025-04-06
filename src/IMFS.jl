module IMFS

using CubicSplines


# difference ẋ[t] = x[t+1] - x[t]
"""
    x[t]  ←  x[t+1]
      ↓         
    ẋ[t]
"""
function bwddif(x::Vector{T}) where T<:AbstractFloat
    N = length(x)
    ẋ = Vector{T}(undef, N)
    ẋ[1:N-1] = x[2:N] - x[1:N-1]
    ẋ[N] = ẋ[N-1]
    return ẋ
end


# fwd difference ẋ[t+1] = x[t+1] - x[t]
"""
    x[t]  →  x[t+1]
                ↓
             ẋ[t+1]
"""
function fwddif(x::Vector{T}) where T<:AbstractFloat
    N = length(x)
    ẋ = Vector{T}(undef, N)
    ẋ[2:N] = x[2:N] - x[1:N-1]
    ẋ[1] = ẋ[2]
    return ẋ
end


"""
    minmaxs(x::Vector{T}) -> y::Vector{T}
return local minimums and local maximums.
    if y[i]<0 , then x[i] is a local maximum.
    if y[i]>0 , then x[i] is a local minimum.
"""
function minmaxs(x::Vector{T}) where T<:AbstractFloat
    return bwddif(sign.(fwddif(x)))
end


export ulenv
"""
    ulenv(x::Vector{T}) -> upper_envelope::Vector{T}, lower_envelope::Vector{T}
"""
function ulenv(x::Vector{T}) where T<:AbstractFloat
    N = length(x)
    t = collect(1:N) .* one(T)
    mms = minmaxs(x)

    # --------------- get upper envelope ids
    mms[1] = -2
    mms[N] = -2
    IMAX = mms .< 0
    # --------------- get lower envelope ids
    mms[1] = 2
    mms[N] = 2
    IMIN = mms .> 0

    U = CubicSpline(t[IMAX], x[IMAX]) # upper envelope
    L = CubicSpline(t[IMIN], x[IMIN]) # lower envelope
    return U[t], L[t]
end


export imfinfos
function imfinfos(x::Vector{T}) where T<:AbstractFloat
    N = length(x)
    t = collect(1:N) .* one(T)
    mms = minmaxs(x)
    ps  = N - sum(iszero.(mms))          # number of extremas
    zc  = sum((x[1:N-1] .* x[2:N]) .< 0) # zero crossing counts
    isimf = abs(ps - zc) ≤ 1
    # --------------- get upper envelope ids
    mms[1] = -2
    mms[N] = -2
    IMAX = mms .< 0
    # --------------- get lower envelope ids
    mms[1] = 2
    mms[N] = 2
    IMIN = mms .> 0

    U = CubicSpline(t[IMAX], x[IMAX]) # upper envelope
    L = CubicSpline(t[IMIN], x[IMIN]) # lower envelope
    M = T(1/2) * (U[t] + L[t])        # middle envelope
    return M, isimf
end


function dB(x::Vector{T}, n::Vector{T}) where T<:AbstractFloat
    ϵ  = T(1e-17)
    ex = sum(x .^ 2) + ϵ
    en = sum(n .^ 2) + ϵ
    return 10log10(ex / en)
end


export nemd
# return n emd
function nemd(x::Vector{T}, N::Int; minsnr::Real=0, maxiters::Int=5) where T<:AbstractFloat
    y = Vector{Vector{T}}(undef,N+1)
    for n = 1:N
        isimf = false
        snr = -Inf
        C = 0
        s = x
        m = x
        while (!isimf) && (snr < minsnr)
            m, isimf = imfinfos(s) # trend and strict imf stop criterion
            r = s - m              # residual = signal - trend
            snr = dB(s, m)         # approx signal to noise ratio in dB
            s = r                  # store residual to s
            # another stop criterion
            C = C + 1
            C ≥ maxiters && break
        end
        y[n] = s    # store residual or so called imf component
        x = x - s   # keep trend
    end
    y[N+1] = x
    return y
end



# 局部中线=上下包络平均，即局部趋势线，抗噪, 当时
# 三次拟合至少需要5个点，所以有可能拟合失败
export localmid
function localmid(x::Vector{T}) where T<:AbstractFloat
    N = length(x)
    t = collect(1:N) .* one(T)
    mms = minmaxs(x)

    # --------------- get upper envelope ids
    mms[1] = -2
    mms[N] = -2
    IMAX = mms .< 0
    # --------------- get lower envelope ids
    mms[1] = 2
    mms[N] = 2
    IMIN = mms .> 0

    U = CubicSpline(t[IMAX], x[IMAX]) # upper envelope
    L = CubicSpline(t[IMIN], x[IMIN]) # lower envelope
    M = T(1/2) * (U[t] + L[t])        # middle envelope
    return M
end


# simple emd
export semd
function semd(x::Vector{T}, N::Int) where T<:AbstractFloat
    y = Vector{Vector{T}}(undef, N+1)
    for n = 1:N
        m = localmid(x)  # trend signal
        r = x - m        # residual signal
        x = m            # trend as new input
        y[n] = r         # store residual
    end
    y[N+1] = x
    return y
end



function locations(x::Vector{T}) where T<:AbstractFloat
    N   = length(x)
    mms = minmaxs(x)
    ps  = N - sum(iszero.(mms))          # number of extremas
    zc  = sum((x[1:N-1] .* x[2:N]) .< 0) # zero crossing counts

    # --------------- get upper envelope ids
    mms[1] = -2
    mms[N] = -2
    IMAX = mms .< 0
    # --------------- get lower envelope ids
    mms[1] = 2
    mms[N] = 2
    IMIN = mms .> 0

    return IMAX, IMIN, abs(ps - zc) ≤ 1
end


function getmid(x::Vector{T}, IMAX::BitVector, IMIN::BitVector) where T<:AbstractFloat
    N = length(x)
    t = collect(1:N) .* one(T)
    U = CubicSpline(t[IMAX], x[IMAX]) # upper envelope
    L = CubicSpline(t[IMIN], x[IMIN]) # lower envelope
    M = T(1/2) * (U[t] + L[t])        # middle envelope
    return M
end


export emd
function emd(x::Vector{T}, N::Int; minsnr::Real=0, maxiters::Int=5) where T<:AbstractFloat
    y = Vector{Vector{T}}()
    for n = 1:N
        isimf = false
        snr = -Inf
        C = 0
        s = x
        m = x
        while (!isimf) && (snr < minsnr)
            imax, imin, isimf = locations(s)
            if sum(imax) < 5 || sum(imin) < 5
                push!(y, x)
                return y
            end
            m = getmid(s, imax,imin) # trend and strict imf stop criterion
            r = s - m                # residual = signal - trend
            snr = dB(s, m)           # approx signal to noise ratio in dB
            s = r                    # store residual to s
            # another stop criterion
            C = C + 1
            C ≥ maxiters && break
        end
        push!(y, s)    # store residual or so called imf component
        x = x - s      # keep trend
    end
    push!(y, x)
    return y
end



end # module IMFS
