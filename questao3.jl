using Plots

include("./metodos.jl")
import .Metodos
include("./equacoes.jl")
import .Burden, .Wikipedia, .Chapra, .Questao3

PATH = "questao3/"

function h_plot(equation, h, file, N = nothing)
    t = Array{Float64}(undef, 0)
    if N === nothing
        N = 0 # qtd de passos ate chegar no final do intervalo
        #t = zeros(Float64, 0)
        for i in equation.ta:h:equation.tb
            N += 1
            append!(t, i)
        end
    else
        for i in 1:N
            append!(t, equation.ta + h*(i-1))
        end
    end

    write(file, "h = $h\nt = $t\n\n")

    y = Metodos.analytic_solution(equation.ta, h, N, t, equation.fy)
    write(file, "Solucao Analitica = $y\n")
    func = plot(t, y, marker=:circle, label="S. Analítica")

    w1_euler, w2_euler = Metodos.euler_sistema(h, N, t, equation.y1_0, equation.y2_0, equation.f1, equation.f2)
    write(file, "Euler: \n\tw1 = $w1_euler\n\tw2 = $w2_euler\n\n")
    plot!(t, w1_euler, marker=:circle, label="Euler")

    w1_rk_4ordem, w2_rk_4ordem = Metodos.rk_4ordem_sistema(h, N, t, equation.y1_0, equation.y2_0, equation.f1, equation.f2)
    write(file, "Runge-Kutta 4 ordem: \n\tw1 = $w1_rk_4ordem\n\tw2 = $w2_rk_4ordem\n\n")
    plot!(t, w1_rk_4ordem, marker=:circle, label="RK 4 ordem")

    w1_adams_moulton, w2_adams_moulton = Metodos.adams_moulton_4_sistema(h, N, t, equation.y1_0, equation.y2_0, equation.f1, equation.f2)
    write(file, "Adams-Moulton 3 passos: \n\tw1 = $w1_adams_moulton\n\tw2 = $w2_adams_moulton\n\n")
    plot!(t, w1_adams_moulton, marker=:circle, label="Adams M. 3 passos")

    filename = string(PATH, "Questao3", "_h_", h, ".png")
    savefig(func, filename)
    
end

function main()
    file = open(string(PATH, "tabela.txt"), "w")

    write(file, "Questao3\n")
    h_plot(Questao3, 1.0, file, 10)
    h_plot(Questao3, 0.5, file, 20)
    h_plot(Questao3, 0.25, file, 40)
    h_plot(Questao3, 0.1, file, 80)

    close(file)
end

main()