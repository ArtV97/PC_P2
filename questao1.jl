using Plots

include("./metodos.jl")
import .Metodos
include("./equacoes.jl")
import .Burden, .Wikipedia, .Chapra


function h_plot(equation, h)
    N = 0 # qtd de passos ate chegar no final do intervalo
    t = zeros(Float64, 0)
    for i in equation.ta:h:equation.tb
        N += 1
        append!(t, i)
    end

    y = Metodos.analytic_solution(equation.ta, h, N, t, equation.fy)
    func = plot(t, y, marker=:circle, label="S. Analítica")

    a2 = 0.5
    w_heun = Metodos.rk_2ordem(a2, h, N, t, equation.y0, equation.f)
    plot!(t, w_heun, marker=:circle, label="Heun")

    a2 = 1
    w_ponto_medio = Metodos.rk_2ordem(a2, h, N, t, equation.y0, equation.f)
    plot!(t, w_ponto_medio, marker=:circle, label="P. Médio")

    a2 = 2/3
    w_ralston = Metodos.rk_2ordem(a2, h, N, t, equation.y0, equation.f)
    plot!(t, w_ralston, marker=:circle, label="Ralston")

    w_euler = Metodos.euler(h, N, t, equation.y0, equation.f)
    plot!(t, w_euler, marker=:circle, label="Euler")
    
    w_rk4 = Metodos.rk_4ordem(h, N, t, equation.y0, equation.f)
    plot!(t, w_rk4, marker=:circle, label="rk 4 ordem")

    filename = string(equation.name, "_h_", h, ".png")
    savefig(func, filename)
end

# h = 2.0/1.0/0.5/0.1
# comparar com as solucoes analiticas
# comparar com Euler explicito e rk_4ordem
function main()
    # Burden
    h_plot(Burden, 2.0)

    h_plot(Burden, 1.0)
    h_plot(Burden, 0.5)
    h_plot(Burden, 0.1)

    # Wikipedia
    h_plot(Wikipedia, 2.0)
    h_plot(Wikipedia, 1.0)
    h_plot(Wikipedia, 0.5)
    h_plot(Wikipedia, 0.1)

    # Chapra
    h_plot(Chapra, 2.0)
    h_plot(Chapra, 1.0)
    h_plot(Chapra, 0.5)
    h_plot(Chapra, 0.1)

end

main()