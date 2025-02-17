function Base.bitreverse(x::Int, n::Int)
    return parse(Int, reverse(bitstring(x)[end - n + 1:end]), base = 2)
end

bit(s::Int, k::Int) = (s >> (k-1)) & 1
flip(s::Int, k::Union{Int,Vector}) = s ⊻ sum(@. 1 << (k-1))

function count1s(i::Int)
    i = i - ((i >>> 1) & 0x55555555)
    i = (i & 0x33333333) + ((i >>> 2) & 0x33333333)
    i = (i + (i >>> 4)) & 0x0f0f0f0f
    i = i + (i >>> 8)
    i = i + (i >>> 16)
    return i & 0x3f
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

function state2bit(index::Int64,N::Int64 = ceil(Int, log2(index))+1)
    str = bitstring(index)[end-N+1:end]
    println(str)
    return str
end

function state2bit(indexs::Union{AbstractVector, AbstractRange},N::Int64 = ceil(Int, log2(indexs[1])) + 1)
    for i in indexs
        println(bitstring(i)[end-N+1:end])
    end
end


function bitpermute(x::Integer, perm::AbstractVector{<:Integer}, N::Integer)
    @assert length(perm) == N "Permutation length must match N"
    @assert sort(perm) == 1:N "Invalid permutation"

    bits_foward = reverse(digits(x, base=2, pad=N))

    new_bits = zeros(Int, N)
    for i in 1:N
        new_bits[perm[i]] = bits_foward[i]
    end

    return foldl((acc, b) -> acc * 2 + b, new_bits; init=0)
end

function statepermute(x::Integer, perm::AbstractVector{<:Integer}, N::Integer, intr::Int64=2)
    @assert length(perm) == N
    state = 0
    for i in 1:intr
        state += bitpermute((x & ((2^N-1) << ((i-1)*N))) >> ((i-1)*N), perm, N) << ((i-1)*N)
    end
    return state
end
function bit2state(str::String)
    return parse(Int64, str, base=2)
end


