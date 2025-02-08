using LinearAlgebra,CairoMakie

function apply_annihilate(j, s)
    if (s & (1 << j)) == 0
        return (0, 0)
    end
    s_new = s ⊻ (1 << j)
    mask = (1 << j) - 1
    num_before = count_ones(s & mask)
    sign = (-1)^num_before
    return (s_new, sign)
end

function apply_create(j, s)
    if (s & (1 << j)) != 0
        return (0, 0)
    end
    s_new = s | (1 << j)
    mask = (1 << j) - 1
    num_before = count_ones(s & mask)
    sign = (-1)^num_before
    return (s_new, sign)
end

function build_hamiltonian(N; t,μ)
    dim = 2^N
    H = zeros(Float64, dim, dim)
    for s in 0:dim-1 # Open boundary condition
        for i in 0:N-2 
            # Term: -t * c_i† c_{i+1}
            (s1, sign1) = apply_annihilate(i+1, s)
            (s2, sign2) = apply_create(i, s1)
            H[s2+1, s+1] += -t * sign1 * sign2

            # Hermitian conjugate term: -t * c_{i+1}† c_i
            (s1_herm, sign1_herm) = apply_annihilate(i, s)
            (s2_herm, sign2_herm) = apply_create(i+1, s1_herm)
            H[s2_herm+1, s+1] += -t * sign1_herm * sign2_herm
        end

        H[s+1,s+1] += -μ*count_ones(s)
    end
    return H
end

function ϵ(k;t,μ)
    return -2*t*cos(k)-μ
end

function getGStheo(N;t,μ)
    k = pi*(1:N) / (N+1)
    ϵs = filter!(x -> x<0, ϵ.(k;t=t,μ = μ))
    E = sum(ϵs)
    N = length(ϵs)
    return E,N
end


# Example usage
N = 4
t = 1.0

lsμ = range(-2,2,51)

Nop = getNop(N)

Es = []
Est = []
Nst = []
Ns = []
for μ in lsμ
    params = (
        t = t,
        μ = μ
    )
    H = build_hamiltonian(N;params...)
    F = eigen(H)
    eigenvalues = F.values
    eigenvectors = F.vectors
    n = eigenvectors[:,1]'*Nop*eigenvectors[:,1]
    push!(Es,eigenvalues[1])
    push!(Ns,n)
    Et,Nt = getGStheo(N;params...)
    if μ == 0
        @show Es[end]
    end
    push!(Est,Et)
    push!(Nst,Nt)
end 

fig = Figure()
ax = Axis(fig[1,1];
title = "Spinless free fermion",
width = 400,
height = 200,
xlabel = L"\mu",
ylabel = L"n")

lines!(lsμ,Nst / N;color = :red,label = "Theory")
scatter!(lsμ,Ns / N;label = "DS-R1 code")
axislegend(ax,position = :lt)
resize_to_layout!(fig)
display(fig)

save("FermiHubbard/figures/spinless_free.png",fig)
