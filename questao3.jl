using Plots

include("./metodos.jl")
import .Metodos
include("./equacoes.jl")
import .Burden, .Wikipedia, .Chapra, .Questao3

PATH = "questao3/"

function h_plot(metodo_nome, metodo::Function, equation, h, file, N = nothing)
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

    w1, w2 = metodo(h, N, t, equation.y1_0, equation.y2_0, equation.f1, equation.f2)
    write(file, "$metodo_nome: \n\tw1 = $w1\n\tw2 = $w2\n\n")
    plot!(t, w1, marker=:circle, label=metodo_nome)

    filename = string(PATH, metodo_nome, "_h_", h, ".png")
    savefig(func, filename)
    
end

function main()
    file = open(string(PATH, "tabela.txt"), "w")

    write(file, "Questao3\n")
    h_plot("Euler", Metodos.euler_sistema, Questao3, 0.1, file, 10)
    h_plot("R.K. 4 ordem", Metodos.rk_4ordem_sistema, Questao3, 0.1, file, 10)
    h_plot("Adams M. 3 passos", Metodos.adams_moulton_4_sistema, Questao3, 0.1, file, 10)

    close(file)

end

main()