function getNop(N)
    Nop = zeros(2^N,2^N)
    for s in 0:2^N-1
        Nop[s+1,s+1] = count_ones(s)
    end
    return Nop
end

