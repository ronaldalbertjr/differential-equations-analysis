### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ ae9b732a-2211-4db4-a28d-13999ed901a7
begin
	using LinearAlgebra
	using Plots
end

# ╔═╡ 94ece1c1-5183-43f3-bdc2-35c27a3cb53d
function GaussianQuadrature(f, x_0, x_f)
    
    w = [1, 1]
    p = [-sqrt(3)/3, sqrt(3)/3]
    
    result = 0
    
    for i in 1:length(p)
        result += w[i]*f((p[i]*(x_f-x_0) + x_f+x_0)/2)*(x_f-x_0)/2
    end
    
    return result
    
end


# ╔═╡ fb6c1c1c-8c2f-4c15-883f-fa83a47ebac4

function GalerkinMethod(f, heat_conduction, h, L, qb, pb)
    x = collect(LinRange(0,L,Int(L/h) + 1))
    
    K = zeros((length(x) - 1, length(x) - 1))
    
    K[1, 1] = (1/h^2)*GaussianQuadrature(heat_conduction, x[1], x[2])
    
    for i in 1:length(x) - 1
        if (i > 1) 
            K[i, i] = ((1/h)^2)*(GaussianQuadrature(heat_conduction, x[i - 1], x[i]) + GaussianQuadrature(heat_conduction, x[i], x[i + 1]))
        end
        if (i < length(x) - 1) 
            K[i, i + 1] = -((1/h)^2)*GaussianQuadrature(heat_conduction, x[i], x[i + 1])
            K[i + 1, i] = K[i, i + 1]
        end
    end
    
    vetor_f = f.(x)
    
    
    F = zeros(length(x) - 1, 1)
    F[1] = (h/3)*vetor_f[1] + (h/6)*vetor_f[2] + qb
    
    for i in 2:length(x) - 1
        F[i] = (h/6)*vetor_f[i - 1] + (2*h/3)*vetor_f[i] + (h/6)*vetor_f[i+1]
    end
    
    F[length(x) - 1] = (h/6)*vetor_f[length(x) - 2] + (2*h/3)*vetor_f[length(x) - 1] + (h/6)*vetor_f[length(x)] + (pb/h^2)*GaussianQuadrature(heat_conduction, x[length(x)-1], x[length(x)])
    
    coef = K\F
    coef = [coef; pb]
    
    return coef

end

# ╔═╡ 1a607e1c-c3bd-4aef-b78b-0676d9706368
begin
# Definindo as constantes do problema
k = 35 
U = 67967200 
L = 0.025 
T_inf = 293.15 
q = 0
h =  0.0025
	

solucao_analitica = x -> T_inf + (U*L^2/(2*k))*(1 - (x/L)^2)
end

# ╔═╡ 55bf85f2-9e34-4140-90d0-9bfa66a08da3
# Gerando os pontos do nosso espaço discretizado
x = collect(LinRange(0,L,Int(L/h) + 1))

# ╔═╡ a95915bb-175c-40e9-8641-c51186383be1
# Temos o resultado da temperatura em cada ponto para cada ponto no vetor x
results = GalerkinMethod(x -> U, x -> k, h, L, q, T_inf)

# ╔═╡ 7e09b891-7859-4723-9806-681443654eb5
# Calculando os resultados reais do vetor x
results_analitica = solucao_analitica.(x)

# ╔═╡ 076b85aa-7f69-4b0f-adab-9451e3760a2b
# Podemos ver que o erro quadrático médio é quase zero.
norm(results_analitica - results)/length(results)

# ╔═╡ 57e9bd22-ec7a-4b02-9ba4-c3692c6b8279
# Os gráficos da solução análitica e dos resultados gerados pelo método numérico, acabaram sobrepostos.
begin
	plot(x,results_analitica,label = "Resultados - Solução Analítica")
	scatter!(x, results, label="Resultados - Método Numérico")
end

# ╔═╡ 30585144-165c-40d0-916e-814980bb7739


# ╔═╡ Cell order:
# ╠═ae9b732a-2211-4db4-a28d-13999ed901a7
# ╠═94ece1c1-5183-43f3-bdc2-35c27a3cb53d
# ╠═fb6c1c1c-8c2f-4c15-883f-fa83a47ebac4
# ╠═1a607e1c-c3bd-4aef-b78b-0676d9706368
# ╠═55bf85f2-9e34-4140-90d0-9bfa66a08da3
# ╠═a95915bb-175c-40e9-8641-c51186383be1
# ╠═7e09b891-7859-4723-9806-681443654eb5
# ╠═076b85aa-7f69-4b0f-adab-9451e3760a2b
# ╠═57e9bd22-ec7a-4b02-9ba4-c3692c6b8279
# ╠═30585144-165c-40d0-916e-814980bb7739
