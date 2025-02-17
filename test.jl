using LatticeUtilities,CairoMakie
include("src/iED.jl")

L = 3
W = 3
N = L*W

Latt = PeriSqua(L,W)

fig = Figure()
ax = Axis(fig[1,1])
latticescatter!(ax,Latt)
display(fig)

displacement(Latt,1,2)
@time dntb = build_destination_array(Latt)
#trantb[1,2,1]
@time nbs,maps = build_neighbor_array(Latt;ordered = true,level=2)
states = collect(0:2^(2*Latt.lattice.N)-1)
mstates = let 
    tmp = []
    for s in states
        CalcMagmom(s,N) == 0 && push!(tmp,s)
    end
    tmp
end

@time transvectb = build_translation_vec_map(Latt,mstates)
#sizeof(transvectb)
tmp = filter!(y -> !isequal(transvectb[y],[3,3]),filter!(x -> length(x) != 2, collect(keys(transvectb))))
vs = map(z -> transvectb[z],tmp)

