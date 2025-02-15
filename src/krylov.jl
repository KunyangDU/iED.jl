
function KrylovTruncate(
                        H::AbstractMatrix,
                        β::Number,
                        Es::Vector,
                        Us::Vector
                        ;
                        max_size = 50,
                        E_tol = 1e-5,
                        projection_scale = 1e3,
                        eigsolve_conditions = (issymmetric=true,ishermitian=true),
                        batch_size = 5,
                        kwargs...
    )

    Es,Us,ϵ,iter = let 

        ϵ = Inf
        iter = length(Es)
        tmpH = deepcopy(H) + projection_scale * (Us[end] * Us[end]')

        while ϵ > E_tol && iter < max_size
            
            F = eigsolve(tmpH,batch_size,:SR;eigsolve_conditions...)

            Es = vcat(Es,F[1][1:batch_size])
            Us = vcat(Us,F[2][1:batch_size])
            iter += batch_size

            ϵ = exp(-β * (Es[end] - Es[1]))
            tmpH .+= sum([projection_scale * (F[2][i] * F[2][i]') for i in 1:batch_size])
            #@show iter,length(Es),batch_size,length(F[1][1:batch_size])

            @assert iter == length(Es)

        end
        Es,Us,ϵ,iter
    end

    return Es,Us,ϵ,iter
end

function KrylovTruncate(H::AbstractMatrix,β::Number,E::Number,U::Number;kwargs...)
    return KrylovTruncate(H,β,[E,],[U,];kwargs...)    
end

function KrylovTruncate(H::AbstractMatrix,β::Number;eigsolve_conditions = (issymmetric=true,ishermitian=true),kwargs...)
    F = eigsolve(H,1,:SR;eigsolve_conditions...)
    return KrylovTruncate(H,β,F[1],F[2];eigsolve_conditions=eigsolve_conditions,kwargs...)
end

function calcEg(H::AbstractMatrix;eigsolve_conditions = (issymmetric=true,ishermitian=true))
    F = eigsolve(H,1,:SR;eigsolve_conditions...)
    E = F[1][1]
    U = F[2][1]
    return E,U
end

