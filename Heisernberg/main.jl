

include("utils.jl")

function main_U1(N = 4)
    state = U1Symm(N, 0)
    H = zeros(2^N,2^N)
    for (ia,a) in enumerate(state)
        for i in 1:N
            j = mod(i,N) + 1
            if bit(a,i) == bit(a,j)
                H[ia,ia] += 1/4
            else
                H[ia,ia] += - 1/4
                b = flip(a,[i,j])
                ib = findfirst(x -> x == b, state)
                H[ia,ib] += 1/2
            end
        end
    end

    return - eigvals(H)[1] / N
end

function main_U1pLatt(N = 4)
    lsk = (-N/2+1):N/2
    state = U1Symm(N,0)
    E0s = Vector(undef,N)
    for (ik,k) in enumerate(lsk)
        pstate,Rs = LattSymm(N, k, state)
        isempty(pstate) && continue
        L = length(pstate)
        H = zeros(L,L) * 1im
        for (ia,a) in enumerate(pstate)
            for i in 1:N
                j = mod(i,N) + 1
                if bit(a,i) == bit(a,j)
                    H[ia,ia] += 1/4
                else
                    H[ia,ia] += - 1/4
                    b = flip(a,[i,j])
                    b,l = LattGauge(b, N)
                    ib = findfirst(x -> x == b, pstate)
                    if !isnothing(ib)
                        H[ia,ib] += 1/2 * sqrt(Rs[ia]/Rs[ib]) * exp(1im*k*l * 2pi / N)
                    end
                end
            end
        end
        H = hermitianize(H)
        E0s[ik] = - eigvals(H)[1] / N
    end
    E0s
end

N = 8
#@time E1 = main_U1(N)
#@time E2 = main_U1pLatt(N)
#@show maximum(E2) ≈ E1;
#E2

#= 
reversion symmetry、spin inversion symmetry的哈密顿量
热容、磁化率的热力学量
=#
state = U1Symm(N, 0)
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

function LattRevSymm(N, k, p)
    ps = []
    Rs = []
    ms = []
    for s in state
        m = nothing
        R = CheckStatePeri(s, N ,k)
        if !isnothing(R)
            R,m = CheckStateRev(s, N, R)
        end
        for σ in [-1,1]
            σ == -1 && k ∈ [0,N/2] && continue
            if !isnothing(m)
                if 1 + σ*p*cos(im*k*m*2pi/N) == 0
                    R = nothing
                end
                if σ == -1 && 1 - σ*p*cos(im*k*m*2pi/N) != 0 
                    R = nothing
                end
            end
            if !isnothing(R)
                push!(ps,s)
                push!(Rs,R)
                push!(ms,m)
            end
        end
    end

    ps,Rs,ms
end
function LattGauge(s, N)
    r = s
    l = 0
    for i in 1:N-1
        t = rotate(s, i, N)
        showstate(t,N)
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

k = 0
p = 1

ps,Rs,ms = LattRevSymm(N, k, p)
LattRevGauge(flip(ps[1],[4,5]),N),ps

#= 
check T
check P
find smallest rep with (R,m) to save σ
=#


