"""
find with julia inbuild function searchsorted(), which is slower than bifind at large array (length > 200).
"""
function ssfind(arr, known_value, known_index, target)
    target == known_value && return known_index
    
    if target > known_value
        range_start = known_index +1
        range_end = lastindex(arr)
    else
        range_start = firstindex(arr)
        range_end = known_index -1
    end

    pos = searchsorted(arr[range_start:range_end], target)

    return length(pos) > 0 ? first(pos) + range_start - 1 : nothing
end
"""
find with binary method, which is faster than ssfind at large array (length > 200).
"""
function bifind(arr, v, k, target)
    if target == v
        return k
    elseif target < v
        left = 1
        right = k - 1
    else
        left = k + 1
        right = length(arr)
    end

    while left <= right
        mid = div((left + right),2)
        if arr[mid] == target
            return mid
        elseif arr[mid] < target
            left = mid + 1
        else
            right = mid - 1
        end
    end
    return nothing
end


function hermitianize(M::Matrix;tol = 1e-5,check=false)
    if check 
        @assert ishermitian(M;tol=tol)
    end
    return M
end

function ishermitian(M::Matrix;tol = 1e-5)
    ϵ = sum(abs.(M .- M'))
    if ϵ < tol
        return true
    else
        return false
    end
end





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

function showstate(index::Int64,N::Int64 = ceil(Int, log2(index))+1)
    println(bitstring(index)[end-N+1:end])
end

function showstate(indexs::Union{AbstractVector, AbstractRange},N::Int64 = ceil(Int, log2(indexs[1])) + 1)
    for i in indexs
        println(bitstring(i)[end-N+1:end])
    end
end
