function latticescatter!(ax::Axis, Latt::SimpleLattice)
    for i in 1:N
        posi = coordinate(Latt,i)
        scatter!(ax,posi...;
        color = :black)
        text!(ax,1.05*posi...;text = "$i")
    end
end