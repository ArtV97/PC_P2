using Plots

include("./metodos.jl")
import .Metodos
include("./equacoes.jl")
import .Burden, .Wikipedia, .Chapra

PATH = "questao2/"

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

    w_adams_bashforth = Metodos.adams_bashforth_4(h, N, t, equation.y0, equation.f)
    write(file, "Adams-Bashforth 4 passos = $w_adams_bashforth\n")
    plot!(t, w_adams_bashforth, marker=:circle, label="Adams B. 4 passos")

    w_adams_moulton = Metodos.adams_moulton_4(h, N, t, equation.y0, equation.f)
    write(file, "Adams-Moulton 4 passos = $w_adams_moulton\n\n")
    plot!(t, w_adams_moulton, marker=:circle, label="Adams M. 4 passos")

    filename = string(PATH, equation.name, "_h_", h, ".png")
    savefig(func, filename)
end

function main()
    file = open(string(PATH, "tabela.txt"), "w")

    # Burden
    write(file, "Burden\n")
    h_plot(Burden, 2.0, file, 6)
    h_plot(Burden, 1.0, file, 6)
    h_plot(Burden, 0.5, file)
    h_plot(Burden, 0.1, file)

    h_plot(Burden, 4.0, file, 10) # teste passo maior

    # Wikipedia
    write(file, "\nWikipedia\n")
    h_plot(Wikipedia, 2.0, file, 6)
    h_plot(Wikipedia, 1.0, file, 6)
    h_plot(Wikipedia, 0.5, file)
    h_plot(Wikipedia, 0.1, file)

    h_plot(Wikipedia, 4.0, file, 10) # teste passo maior

    # Chapra
    write(file, "\nChapra\n")
    h_plot(Chapra, 2.0, file, 6)
    h_plot(Chapra, 1.0, file, 6)
    h_plot(Chapra, 0.5, file)
    h_plot(Chapra, 0.1, file)

    h_plot(Chapra, 12000.0, file, 10) # teste passo maior

    close(file)

end

main()