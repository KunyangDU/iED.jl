
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




