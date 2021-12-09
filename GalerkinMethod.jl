using LinearAlgebra

function GaussianQuadrature(f, x_0, x_f)
    w = [1, 1]
    p = [-sqrt(3)/3, sqrt(3)/3]
    
    result = 0
    
    for i in 1:length(p)
        result += w[i]*f((p[i]*(x_f-x_0) + x_f+x_0)/2)*(x_f-x_0)/2
    end
    
    return result
end


function GalerkinMethod(f, heat_conduction, h, L, qb, pb)
    x = collect(LinRange(0,L,Int(L/h) + 1))
    
    K = zeros((length(x) - 1, length(x) - 1))
    
    K[1, 1] = (1/h)*GaussianQuadrature(heat_conduction, x[1], x[2])
    
    for i in 1:length(x) - 1
        if (i > 1) 
            K[i, i] = ((1/h)^2)*(GaussianQuadrature(heat_conduction, x[i - 1], x[i]) + GaussianQuadrature(heat_conduction, x[i], x[i + 1]))
        end
        if (i < length(x) - 1) 
            K[i, i + 1] = (-(1/h)^2)*GaussianQuadrature(heat_conduction, x[i], x[i + 1])
            K[i + 1, i] = K[i, i + 1]
        end
    end
    
    vetor_f = f.(x)
    
    
    F = zeros(length(x) - 1, 1)
    F[1] = (h/3)*vetor_f[1] + (h/6)*vetor_f[2] + qb
    
    for i in 2:length(x) - 2
        F[i] = (h/6)*vetor_f[i - 1] + (2*h/3)*vetor_f[i] + (h/6)*vetor_f[i+1]
    end
    
    F[length(x) - 1] = (h/6)*vetor_f[length(x) - 2] + (2*h/3)*vetor_f[length(x) - 1] + (h/6)*vetor_f[length(x)] + (pb/h^2)*GaussianQuadrature(heat_conduction, x[length(x)-1], x[length(x)])
    
    return F\K

end

GalerkinMethod(x -> x^2, x -> x*(1-x) + 1, 0.2, 1, 0, -1)