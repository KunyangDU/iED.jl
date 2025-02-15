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


