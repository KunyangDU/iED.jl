function ϵ(k)
    return -2sum(cos.(k))
end

function getk(L::Int;condition = :obc)
    if condition == :obc
        return @. pi * (1:L) / (L+1)
    elseif condition == :pbc
        return @. 2pi * (1:L) / L
    end
end

function getk(Lx::Int,Ly::Int)
    if Ly == 1
        lsk = getk(Lx)
    else
        lskx = getk(Lx)
        lsky = getk(Ly;condition = :pbc)
        lsk = [[kx,ky] for kx in lskx,ky in lsky][:]
    end
    return lsk
end

function ue(β::Number,Lx::Int,Ly::Int)
    lsk = getk(Lx,Ly)
    lsum = @.  ϵ(lsk) / (1 + exp( β * ϵ(lsk)))
    return sum(lsum) / Lx / Ly
end

function fe(β::Number,Lx::Int,Ly::Int)
    lsk = getk(Lx,Ly)
    return - sum(@. log(1+exp(-β*(ϵ(lsk))))) / β / Lx / Ly
end

function ce(β::Number,Lx::Int,Ly::Int)
    lsk = getk(Lx,Ly)
    return β^2/2 * sum(@. ϵ(lsk)^2/(1 + cosh(β * ϵ(lsk)))) / Lx / Ly
end