using LinearAlgebra,CairoMakie

include("utils.jl")
include("../Heisenberg/utils.jl")

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

function build_hamiltonian(N,states = collect(0:2^N -1); t,μ)
    H = zeros(Float64, length(states), length(states))
    for (ind,s) in enumerate(states) # Open boundary condition
        for i in 0:N-2
            # Term: -t * c_i† c_{i+1}
            (s1, sign1) = apply_annihilate(i+1, s)
            (s2, sign2) = apply_create(i, s1)
            ind2 = findfirst(x -> x == s2,states)

            # Hermitian conjugate term: -t * c_{i+1}† c_i
            (s1_herm, sign1_herm) = apply_annihilate(i, s)
            (s2_herm, sign2_herm) = apply_create(i+1, s1_herm)
            ind2_herm = findfirst(x -> x == s2_herm,states)

            !isnothing(ind2) && (H[ind2, ind] += -t * sign1 * sign2)
            !isnothing(ind2_herm) && (H[ind2_herm, ind] += -t * sign1_herm * sign2_herm)
        end

        H[ind,ind] += -μ*count_ones(s)
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
N = 9
t = 1.0

lsμ = range(-2,2,51)

Nop = getNop(N)

Es = []
Est = []
Nst = []
Ns = []
lsβ = 10 .^ range(log10.([5e-3,50])...,30)
lsT = 1 ./ lsβ

for μ in 0
    params = (
        t = t,
        μ = μ
    )
    data = Dict()

    for N0 in 0:N
        states = let 
            tmp = []
            for s in 0:2^N-1
                count_ones(s) == N0 && push!(tmp,s)
            end
            tmp
        end
        H = build_hamiltonian(N,states;params...)
        F = eigen(H)
        eigenvalues = F.values
        eigenvectors = F.vectors
        tmpdata = Dict(
            "E" => eigenvalues,
            "V" => eigenvectors,
            "N" => N0,
        )
        data[(μ,N0)] = tmpdata
        #n = eigenvectors[:,1]'*Nop*eigenvectors[:,1]
        #push!(Es,eigenvalues[1])
        #push!(Ns,n)
        #Et,Nt = getGStheo(N;params...)
        #@show Es[end]
    end

    Fs = zeros(length(lsβ))
    Us = zeros(length(lsβ))
    Ns = zeros(length(lsβ))
    Ces = zeros(length(lsβ))

    for (iβ,β) in enumerate(lsβ)
        Z = sum(vcat([exp.(- β * data[key]["E"]) for key in keys(data)]...))
        U = sum(vcat([exp.(- β * data[key]["E"]) .* data[key]["E"] for key in keys(data)]...)) / Z
        Ntmp = sum(vcat([exp.(- β * data[key]["E"]) .* data[key]["N"] for key in keys(data)]...)) / Z
        U2 = sum(vcat([exp.(- β * data[key]["E"]) .* (data[key]["E"] .^ 2) for key in keys(data)]...)) / Z

        Fs[iβ] = -log(Z) / β / N
        Us[iβ] = U / N
        Ns[iβ] = Ntmp / N
        Ces[iβ] = (U2 - U^2) * β ^ 2 / N
    end
    figsize = (width = 300,
    height = 200,)

    fig = Figure()
    axf = Axis(fig[1,1];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    scatter!(axf,lsT,Fs)

    axU = Axis(fig[1,2];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    scatter!(axU,lsT,Us)

    axN = Axis(fig[2,1];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    ylims!(axN,0.4,0.6)

    scatter!(axN,lsT,Ns)

    axCe = Axis(fig[2,2];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    scatter!(axCe,lsT,Ces)

    resize_to_layout!(fig)
    display(fig)


end 

