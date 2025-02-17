function neighbor(Latt::SquareLattice;level::Int64 = 1,ordered::Bool = false)
    return collect(Tuple.(sort.(eachcol(build_neighbor_table(Latt.bonds[(ordered,level)], Latt.unitcell, Latt.lattice)))))
end
"""
translation vector on lattice basis
"""
function displacement(Latt::SimpleLattice,i::Int64,j::Int64)
    return collect(sites_to_displacement(i,j,Latt.unitcell,Latt.lattice))
end
"""
translation vector in Descartes coordinate
"""
function distance(Latt::SimpleLattice,i::Int64,j::Int64)
    return collect(displacement_to_vec(displacement(Latt,i,j),1,1,Latt.unitcell))
end
"""
position on lattice basis with sublattice index
"""
function location(Latt::SimpleLattice,i::Int64)
    return site_to_loc(i,Latt.unitcell,Latt.lattice)
end
"""
position in Descartes coordinate
"""
function coordinate(Latt::SimpleLattice,i::Int64)
    return loc_to_pos(location(Latt,i)[1],1,Latt.unitcell)
end
"""
return tb = {L*W*N array}, satisfying tb[i,j,site1] -> site2 which can be obstain by a translation vector [i,j] from site1.
"""
function build_destination_array(Latt::SimpleLattice)
    W,L = Latt.lattice.L
    N = Latt.lattice.N
    tb = zeros(Int64,L,W,N)
    for s in 1:N, i in 1:L,j in 1:W
        tb[i,j,s] = site_to_site(s,[j,i],1,Latt.unitcell,Latt.lattice)
    end
    @assert 0 ∉ tb
    return tb
end
"""
return nbls = {level-th neighbor list}, maps = {array}, satisfying maps[i] = (bonds connecting i, sites level-the neighboring i)
"""
function build_neighbor_array(Latt::SimpleLattice;level::Int64 = 1,ordered::Bool = false)
    maps = map_neighbor_table(build_neighbor_table(Latt.bonds[(ordered,level)],Latt.unitcell,Latt.lattice))
    maps = [maps[i] for i in 1:Latt.lattice.N]
    nbls = neighbor(Latt;level = level,ordered = ordered)
    return nbls,maps
end
"""
return tb = {L*W*length(states) array}, satisfying tb[i,j,state1] -> state2 which can be obstain by a translation vector [i,j] from state1.
"""
function build_translation_array(Latt::SimpleLattice,states::Vector,intr::Int64=2)
    W,L = Latt.lattice.L
    N = Latt.lattice.N
    dntb = build_destination_array(Latt)
    tb = zeros(Int64,L,W,length(states))
    for i in 1:L, j in 1:W 
        perm = dntb[i,j,:]
        for (is,s) in enumerate(states)
            tb[i,j,is] = statepermute(s,perm,N,intr)
        end
    end
    return tb
end
"""
return tb = {Dict}, satisfying:
- tb[(i,)] -> vector which leads i again, i.e., the smallest periodicity; 
- tb[(i,j)] -> vector which leads j , i.e., the translation vector.
"""
function build_translation_vec_map(Latt::SimpleLattice,states::Vector,intr::Int64=2)
    W,L = Latt.lattice.L
    N = Latt.lattice.N
    dntb = build_destination_array(Latt)
    tb = Dict()
    for i in 1:L, j in 1:W
        perm = dntb[i,j,:]
        for s1 in states
            s2 = statepermute(s1,perm,N,intr)
            if s1 == s2
                exi = get(tb,(s1,),nothing)
                tb[(s1,)] = isnothing(exi) ? [i,j] : (sum(exi) < i+j ? exi : [i,j])
            else
                tb[Tuple(sort([s1,s2]))] = [i,j]
            end
        end
    end
    return tb
end

