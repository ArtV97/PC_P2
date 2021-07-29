using Plots

include("./metodos.jl")
import .Metodos
include("./equacoes.jl")
import .Burden, .Wikipedia, .Chapra

PATH = "results/"

function h_plot(equation, h, file)
    N = 0 # qtd de passos ate chegar no final do intervalo
    t = zeros(Float64, 0)
    for i in equation.ta:h:equation.tb
        N += 1
        append!(t, i)
    end

    write(file, "h = $h\nt = $t\n\n")

    y = Metodos.analytic_solution(equation.ta, h, N, t, equation.fy)
    write(file, "Solucao Analitica = $y\n")
    func = plot(t, y, marker=:circle, label="S. Analítica")

    a2 = 0.5
    w_heun = Metodos.rk_2ordem(a2, h, N, t, equation.y0, equation.f)
    write(file, "Heun = $w_heun\n")
    plot!(t, w_heun, marker=:circle, label="Heun")

    a2 = 1
    w_ponto_medio = Metodos.rk_2ordem(a2, h, N, t, equation.y0, equation.f)
    write(file, "Ponto Medio = $w_ponto_medio\n")
    plot!(t, w_ponto_medio, marker=:circle, label="P. Médio")

    a2 = 2/3
    w_ralston = Metodos.rk_2ordem(a2, h, N, t, equation.y0, equation.f)
    write(file, "Ralston = $w_ralston\n")
    plot!(t, w_ralston, marker=:circle, label="Ralston")

    w_euler = Metodos.euler(h, N, t, equation.y0, equation.f)
    write(file, "Euler = $w_euler\n")
    plot!(t, w_euler, marker=:circle, label="Euler")
    
    w_rk4 = Metodos.rk_4ordem(h, N, t, equation.y0, equation.f)
    write(file, "Runge Kutta 4 ordem = $w_rk4\n\n")
    plot!(t, w_rk4, marker=:circle, label="rk 4 ordem")

    filename = string(PATH, equation.name, "_h_", h, ".png")
    savefig(func, filename)
end

# h = 2.0/1.0/0.5/0.1
# comparar com as solucoes analiticas
# comparar com Euler explicito e rk_4ordem
function main()
    file = open(string(PATH, "tabela_questao1.txt"), "w")

    # Burden
    write(file, "Burden\n")
    h_plot(Burden, 2.0, file)
    h_plot(Burden, 1.0, file)
    h_plot(Burden, 0.5, file)
    h_plot(Burden, 0.1, file)

    # Wikipedia
    write(file, "\nWikipedia\n")
    h_plot(Wikipedia, 2.0, file)
    h_plot(Wikipedia, 1.0, file)
    h_plot(Wikipedia, 0.5, file)
    h_plot(Wikipedia, 0.1, file)

    # Chapra
    write(file, "\nChapra\n")
    h_plot(Chapra, 2.0, file)
    h_plot(Chapra, 1.0, file)
    h_plot(Chapra, 0.5, file)
    h_plot(Chapra, 0.1, file)

    close(file)

end

main()