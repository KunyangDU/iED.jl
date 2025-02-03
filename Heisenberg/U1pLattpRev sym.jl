
include("utils.jl")

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

function main_U1pLattpRev(N)
    k = 0
    p = 1

    ps,Rs,ms,σs = LattRevSymm(N, k, p)
    #@show ps
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
            #@show b,l,q
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
                        nm = (1/2) * (σs[i1]*p)^q * (Norm(N,k,Rs[j1],σs[j1],p,ms[j1])/Norm(N,k,Rs[i1],σs[i1],p,ms[i1]))
                        c = let 
                            if isnothing(ms[j1])
                                if σs[i1] == σs[j1]
                                    #@show 1
                                    cos(k*l*2pi/N)
                                else
                                    #@show 2
                                    σs[j1]*sin(k*l*2pi/N)
                                end
                            else
                                #@show ps[j1]
                                #@show l,ms[j1],l-ms[j1]
                                if σs[i1] == σs[j1]
                                    #@show 3
                                    (cos(k*l*2pi/N) + σs[j1]*p*cos(k*(l-ms[j1])*2pi/N)) / (1 + σs[j1]*p*cos(k*ms[j1]*2pi/N))
                                else
                                    σs[j1]*(sin(k*l*2pi/N) + σs[j1]*p*sin(k*(l-ms[j1])*2pi/N)) / (1 + σs[j1]*p*cos(k*ms[j1]*2pi/N))
                                end
                            end    
                        end
                        nm*c
                    end
                    #showstate(a,N)
                    #showstate(b,N)
                    H[i1,j1] += term
                end
            end
        end
    end    
end


#= 
N_b / N_a的计算有误，例如N=8，k=0，p=1的(1,2),(2,1)
=#


