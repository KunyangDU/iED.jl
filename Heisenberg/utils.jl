using LinearAlgebra

mutable struct IndexStates{B, N, L}
    index::Vector 
    mask::Int
    function IndexStates(index::Vector, state::Vector; base =2 , Nmask = ceil(Int, log(base, length(index))))
        @assert length(index) == length(state)
        return new{base, Nmask, length(index)}(index, state, base^N - 1)
    end
end

function rotate(x::Int, n::Int64, Nmask::Int64; direction = :left)
    mask = 2^Nmask - 1
    x = x & mask
    n = n % Nmask
    if direction == :left 
        return ((x << n) | (x >> (Nmask - n))) & mask
    else
        return ((x >> n) | (x << (Nmask - n))) & mask
    end
end

function showstate(index::Int64,N::Int64 = ceil(Int, log2(index))+1)
    println(bitstring(index)[end-N+1:end])
end

function showstate(indexs::Union{AbstractVector, AbstractRange},N::Int64 = ceil(Int, log2(indexs[1])) + 1)
    for i in indexs
        println(bitstring(i)[end-N+1:end])
    end
end

function U1Symm(N, m, state = 0:2^N-1)
    tmp = []
    for s in state
        if (2*count1s(s & (2^N - 1)) - N) // 2 == m 
            push!(tmp, s)
        end
    end
    return tmp
end

function count1s(i::Int)
    i = i - ((i >>> 1) & 0x55555555)
    i = (i & 0x33333333) + ((i >>> 2) & 0x33333333)
    i = (i + (i >>> 4)) & 0x0f0f0f0f
    i = i + (i >>> 8)
    i = i + (i >>> 16)
    return i & 0x3f
end

function periodize(s::Int, k::Number)
    ps = rotate.(s,1:N,N)
    sp = minimum(ps)
    R = unique(diff(findall(x -> x == s,repeat(ps,2))))
    @assert length(R) == 1
    return sp, R[1]
end

bit(s::Int, k::Int) = (s >> (k-1)) & 1
flip(s::Int, k::Union{Int,Vector}) = s ⊻ sum(@. 1 << (k-1))

"""
return periodicity R
"""
function CheckStatePeri(s,N,k)
    for i in 1:N
        t = rotate(s,i,N)
        if t<s
            return nothing 
        elseif t == s 
            if mod(i*k,N) ≠ 0
                return nothing
            end
            return i
        end
    end
end

function LattSymm(N, k, state = 0:2^N-1)
    ps = []
    Rs = []
    for s in state
        R = CheckStatePeri(s,N,k)
        if !isnothing(R)
            push!(ps,s)
            push!(Rs,R)
        end
    end
    return ps,Rs
end

function LattGauge(s, N)
    r = s
    l = 0
    for i in 1:N-1
        t = rotate(s, i, N)
        if t < r
            r = t 
            l=i
        end
    end
    return r,l
end

function hermitianize(M::Matrix;tol = 1e-5)
    ϵ = sum(abs.(M .- M'))
    @assert ϵ < tol
    return (M .+ M') / 2
end

function Base.bitreverse(x::Int, n::Int)
    return parse(Int, reverse(bitstring(x)[end - n + 1:end]), base = 2)
end
"""
return reversion periodicity m
"""
function CheckStateRev(s, N, R)
    t = bitreverse(s,N)
    for i in 0:R-1
        t < s && return nothing, nothing
        t == s && return R, i
        t = rotate(t,1,N)
    end
    return R,nothing
end

function LattGauge(s, N)
    r = s
    l = 0
    for i in 1:N-1
        t = rotate(s, i, N)
        if t < r
            r = t 
            l=i
        end
    end
    return r,l
end


function LattRevGauge(s, N)
    r,l = LattGauge(s,N)
    t = bitreverse(s,N)
    q = 0
    for i in 1:N-1
        t = rotate(t,1,N)
        if t < r 
            r=t 
            l=i 
            q=1
        end
    end 
    return r,l,q
end

function DiagTerm(state::Vector,N::Int)
    ds = zeros(length(state))
    for (is,s) in enumerate(state)
        for i in 1:N
            j = mod(i,N) + 1
            if bit(s,i) == bit(s,j)
                ds[is] += 1/4
            else
                ds[is] += -1/4
            end
        end
    end
    return ds
end

function _getSubsize(ind::Int, state::Vector)
    if ind > 1 && state[ind] == state[ind-1]
        return nothing 
    elseif ind < length(state) && state[ind] == state[ind+1]
        return 2
    else
        return 1
    end
end

