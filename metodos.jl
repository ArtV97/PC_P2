module Metodos

    export analytic_solution, euler, rk_4ordem, rk_2ordem

    function analytic_solution(ta, tb, h, fy::Function)
        N = 0 # qtd de passos ate chegar no final do intervalo
        for i in ta:h:tb
            N += 1
        end
        
        y = zeros(Float64, N)
        t = zeros(Float64, N)

        t[1] = ta
        y[1] = fy(ta)
        # solucao analitica
        for i=1:N
            t[i+1] = a + h*i
            y[i+1] = fy(t[i+1])
        end

        return y
    end

    function euler(ta, tb, h, y0, f::Function)
        N = 0 # qtd de passos ate chegar no final do intervalo
        for i in ta:h:tb
            N += 1
        end
        
        w = zeros(Float64, N)
        t = zeros(Float64, N)

        t[1] = ta
        w[1] = y0
        for i=1:N
            w[i+1] = w[i] + h*f(t[i], w[i])
            t[i+1] = a + h*i
        end

        return w
    end

    function rk_4ordem(ta, tb, h, y0, f::Function)
        N = 0 # qtd de passos ate chegar no final do intervalo
        for i in ta:h:tb
            N += 1
        end

        w = zeros(Float64, N)
        t = zeros(Float64, N)

        t[1] = ta
        w[1] = y0
        for i=1:N
            k1 = f(t[i], w[i])
            k2 = f(t[i]+h/2, w[i]+k1*h/2)
            k3 = f(t[i]+h/2, w[i]+k2*h/2)
            k4 = f(t[i]+h, w[i]+k3*h)
            w[i+1] = w[i] + h*(k1+2*k2+2*k3+k4)/6
        end

        return w
    end

    # Runge-Kutta 2 ordem

    # w0 = alfa
    # wi1 = wi + (a1+k1 + a2k2)*h
    # k1 = f(ti,wi)
    # k2 = f(ti + p1*h, wi + q11*k1*h)

    # a1+a2 = 1
    # a2p1 = 1/2
    # a2q11 = 1/2

    function rk_2ordem(a2, h, ta, tb, y0, f::Function)
        a1 = 1 - a2
        p1 = 0.5 / a2
        q11 = 0.5 / a2

        N = 0 # qtd de passos ate chegar no final do intervalo
        for i in ta:h:tb
            N += 1
        end
        
        w = zeros(Float64, N)
        t = zeros(Float64, N)

        t[1] = ta
        w[1] = y0

        for i in 1:N
            k1 = f(t[i], w[i])
            k2 = f(t[i] * p1*h, wi + q11*k1*h)
            w[i+1] = w[i] + (a1+k1 + a2k2) * h
            
            t[i+1] = a + h*i
        end
        
        return w
    end
end
#=
# rk_2ordem eh usado apenas como "preditor"
function rk_2ordem(a2, h, ti, wi, f::Function)
    a1 = 1 - a2
    p1 = 0.5 / a2
    q11 = 0.5 / a2

    k1 = f(ti, wi)
    k2 = f(ti * p1*h, wi + q11*k1*h)

    return wi + (a1+k1 + a2k2) * h
end

function heun(y0, ta, tb, h, a2, f::Function)
    N = 0 # qtd de passos ate chegar no final do intervalo
    for i in ta:h:tb
        N += 1
    end
    
    w = zeros(Float64, N)
    t = zeros(Float64, N)

    t[1] = ta
    w[1] = y0
    for i in 1:N
        w[i+1] = rk_2ordem(a2, h, t[i], w[i], f)
        t[i+1] = a + h*i
        w[i+1] = w[i] + h*f(t[i+1], w[i+1])
    end
    
    return w
end

function ponto_medio(y0, ta, tb, h, a2, f::Function)
    N = 0 # qtd de passos ate chegar no final do intervalo
    for i in ta:h:tb
        N += 1
    end
    
    w = zeros(Float64, N)
    t = zeros(Float64, N)

    t[1] = ta
    w[1] = y0
    for i in 1:N
        w[i+1] = rk_2ordem(a2, h, t[i], w[i], f)
        t[i+1] = a + h*i

        # inclinacao_media = (derivada1 + derivada2) / 2
        d2 = f(t[i+1], w[i+1])
        derivada_media = (f(t[i], w[i]) + d2) / 2
        w[i+1] = w[i] + h * derivada_media
    end
    
    return w
end
=#