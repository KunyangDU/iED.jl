

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

function OffdiagTerm(i,j,l,q)
    return 0
end

function LattRevSymm(N, k, p)
    ps = []
    Rs = []
    ms = []
    σs = []
    for s in state
        
        m = nothing
        R = CheckStatePeri(s, N ,k)
        if !isnothing(R)
            R,m = CheckStateRev(s, N, R)
        end
        R0 = R
        for σ in [1,-1]
            R = R0
            σ == -1 && k ∈ [0,N/2] && continue
            if !isnothing(m)
                if 1 + σ*p*cos(k*m*2pi/N) == 0
                    R = nothing
                end
                if σ == -1 && 1 - σ*p*cos(k*m*2pi/N) != 0 
                    R = nothing
                end
            end
            if !isnothing(R)
                push!(ps,s)
                push!(Rs,R)
                push!(ms,m)
                push!(σs,σ)
            end
        end
    end

    ps,Rs,ms,σs
end

function Norm(N,k,R,σ,p,m)
    if k ∈ [0,N/2]
        g=1
    else
        g=2
    end
    if isnothing(m)
        return N^2*g/R 
    else
        return N^2*g*(1+σ*p*cos(k*m*2pi/N))/R
    end    
end

k = 1
p = 1

ps,Rs,ms,σs = LattRevSymm(N, k, p)
@show ps
ds = DiagTerm(ps,N)

H = zeros(length(ps),length(ps))

for (ia,a) in enumerate(ps)
    n = _getSubsize(ia,ps)
    isnothing(n) && continue

    for i in ia:ia+n-1
        H[i,i] += ds[i]
    end

    for i in 1:N 
        j = mod(i,N) + 1
        bit(a,i) == bit(a,j) && continue
        t = flip(a,[i,j])
        b,l,q = LattRevGauge(t, N)
        ib = findfirst(x -> x == b, ps)
        @show b,l,q
        if !isnothing(ib)
            if ib > 1 && ps[ib] == ps[ib-1]
                m = 2
                ib -= 1
            elseif ib < length(ps) && ps[ib] == ps[ib+1]
                m = 2
            else
                m = 1
            end

            for i1 in ia:ia+n-1, j1 in ib:ib+m-1
                term = let 
                    nm = (1/2) * (σs[i1]*p)^q * /(map(x -> Norm(N,k,Rs[x],σs[x],p,ms[x]),[j1,i1])...)
                    c = let 
                        if isnothing(ms[j1])
                            if σs[i1] == σs[j1]
                                @show 1
                                cos(k*l*2pi/N)
                            else
                                @show 2
                                σs[j1]*sin(k*l*2pi/N)
                            end
                        else
                            if σs[i1] == σs[j1]
                                @show 3
                                (cos(k*l*2pi/N) + σs[j1]*p*cos(k*(l-ms[j1])*2pi/N)) / (1 + σs[j1]*p*cos(k*ms[j1]*2pi/N))
                            else
                                @show 4
                                σs[j1]*(sin(k*l*2pi/N) + σs[j1]*p*sin(k*(l-ms[j1])*2pi/N)) / (1 + σs[j1]*p*cos(k*ms[j1]*2pi/N))
                            end
                        end    
                    end
                    @show nm,c
                    nm*c
                end
                @show i1,j1,term
                H[i1,j1] += term
            end
        end
    end
end

H
#eigvals(H)[1]

#= 
Hamiltonian for given state ps
=#

#= 
分块化有问题
没有成功按sigma子空间划分
=#


