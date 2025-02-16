using LatticeUtilities,CairoMakie
include("src/iED.jl")
abstract type AbstractLattice end
abstract type SimpleLattice <: AbstractLattice end
mutable struct SquareLattice <: SimpleLattice
    unitcell::LatticeUtilities.UnitCell
    lattice::LatticeUtilities.Lattice
    bonds::Dict
    function SquareLattice(unitcell,lattice,bonds)
        return new(unitcell,lattice,bonds)
    end
end

function PeriSqua(L,W)
    square = UnitCell(lattice_vecs = [[0.,1.],[1.,0.]],basis_vecs = [[0.,0.]])
    lattice = Lattice(L = [W,L], periodic = [true,true])
    bonds = Dict(
        (true,1) => vcat([[Bond((1,1), [i,0]),Bond((1,1), [0,i])] for i in [-1,1]]...),
        (false,1) => [Bond((1,1), [1,0]),Bond((1,1), [0,1])],
        (true,2) => [Bond((1,1), [i,j]) for i in [-1,1],j in [-1,1]][:],
        (false,2) => [Bond((1,1), [1,i]) for i in [-1,1]],
    )
    return SquareLattice(square,lattice,bonds)
end

function latticescatter!(ax::Axis, Latt::SimpleLattice)
    for i in 1:N
        posi = coordinate(Latt,i)
        scatter!(ax,posi...;
        color = :black)
        text!(ax,1.05*posi...;text = "$i")
    end
end

function neighbor(Latt::SquareLattice;level::Int64 = 1,ordered::Bool = false)
    return collect(Tuple.(sort.(eachcol(build_neighbor_table(Latt.bonds[(ordered,level)], Latt.unitcell, Latt.lattice)))))
end

function build_destination_table(Latt::SimpleLattice)
    W,L = Latt.lattice.L
    N = Latt.lattice.N
    tb = zeros(Int64,L,W,N)
    for s in 1:N, i in 1:L,j in 1:W
        tb[i,j,s] = site_to_site(s,[j,i],1,Latt.unitcell,Latt.lattice)
    end
    @assert 0 ∉ tb
    return tb
end

function displacement(Latt::SimpleLattice,i::Int64,j::Int64)
    return collect(sites_to_displacement(i,j,Latt.unitcell,Latt.lattice))
end

function distance(Latt::SimpleLattice,i::Int64,j::Int64)
    return collect(displacement_to_vec(displacement(Latt,i,j),1,1,Latt.unitcell))
end

function coordinate(Latt::SimpleLattice,i::Int64)
    return loc_to_pos(site_to_loc(i,Latt.unitcell,Latt.lattice)[1],1,Latt.unitcell)
end

function build_neighbor_map(Latt::SimpleLattice;level::Int64 = 1,ordered::Bool = false)
    maps = map_neighbor_table(build_neighbor_table(Latt.bonds[(ordered,level)],Latt.unitcell,Latt.lattice))
    maps = [maps[i] for i in 1:Latt.lattice.N]
    nbls = neighbor(Latt;level = level,ordered = ordered)
    return nbls,maps
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
    state = 0
    for i in 1:intr
        state += bitpermute((x & ((2^N-1) << ((i-1)*N))) >> ((i-1)*N),perm,N) << ((i-1)*N)
    end
    return state
end

function build_translation_table(Latt::SimpleLattice,states::Vector,intr::Int64=2)
    W,L = Latt.lattice.L
    N = Latt.lattice.N
    dntb = build_destination_table(Latt)
    tb = zeros(Int64,L,W,length(states))
    for i in 1:L,j in 1:W 
        perm = dntb[i,j,:]
        for (is,s) in enumerate(states)
            tb[i,j,is] = statepermute(s,perm,N,intr)
        end
    end
    return tb
end


L = 4
W = 4
N = L*W

Latt = PeriSqua(L,W)

fig = Figure()
ax = Axis(fig[1,1])
latticescatter!(ax,Latt)
display(fig)

displacement(Latt,1,2)
dntb = build_destination_table(Latt)
#trantb[1,2,1]
nbs,maps = build_neighbor_map(Latt;ordered = true,level=2)
dntb[1,2,:]

build_translation_table(Latt,collect(0:2^Latt.lattice.N-1))