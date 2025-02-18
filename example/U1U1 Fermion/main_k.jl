include("../../src/iED.jl")

function U1filter(filt::Function,Latt::SimpleLattice,u::Number,states::Vector = collect(0:2^(2*Latt.lattice.N)-1))
    u1states = let 
        tmp = []
        for s in states
            filt(s,N) == u && push!(tmp,s)
        end
        tmp
    end
    return u1states
end

function CalcPartnum(s::Int64,::Int64)
    return count1s(s)
end

L = 4
W = 4
N = L*W

Latt = PeriSqua(L,W)

nstates = U1filter(CalcPartnum,Latt,N)
nmstates = U1filter(CalcMagmom,Latt,0,nstates)


