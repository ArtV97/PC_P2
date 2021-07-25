using Plots

include("./metodos.jl")
import .Metodos
include("./equacoes.jl")
import .Burden, .Wikipedia, .Chapra

#Metodos.teste()
#println(Burden.f(0,1))
#println("ta: ", Burden.ta)
#println("tb: ", Burden.tb)
println(Wikipedia.fy(3))

# h = 2.0/1.0/0.5/0.1
# comparar com as solucoes analiticas
# comparar com Euler explicito e rk_4ordem

# Burden
y = Metodos.analytic_solution(Burden.ta, Burden.tb, Burden.h, Burden.fy)


