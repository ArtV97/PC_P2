module Metodos

    export analytic_solution, euler, rk_4ordem, rk_2ordem, adams_bashforth_4, adams_moulton_4, euler_sistema, rk_4ordem_sistema, adams_moulton_4_sistema

    function analytic_solution(ta, h, N, t::Array, fy::Function)        
        y = zeros(Float64, N)
        y[1] = fy(ta)
        # solucao analitica
        for i=1:N-1
            y[i+1] = fy(t[i+1])
        end

        return y
    end

    function euler(h, N, t, y0, f::Function)        
        w = zeros(Float64, N)
        w[1] = y0

        for i=1:N-1
            w[i+1] = w[i] + h*f(t[i], w[i])
        end

        return w
    end

    function rk_4ordem(h, N, t, y0, f::Function)
        w = zeros(Float64, N)
        w[1] = y0

        for i=1:N-1
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

    function rk_2ordem(a2, h, N, t, y0, f::Function)
        a1 = 1 - a2
        p1 = 0.5 / a2
        q11 = 0.5 / a2
        
        w = zeros(Float64, N)
        w[1] = y0

        for i in 1:N-1
            k1 = f(t[i], w[i])
            k2 = f(t[i] + p1*h, w[i] + q11*k1*h)
            w[i+1] = w[i] + (a1*k1 + a2*k2) * h
        end
        
        return w
    end

    # Adams-Bashforth 4 ordem 4 passos
    function adams_bashforth_4(h, N, t, y0, f::Function)
        w = zeros(Float64, N)
        w[1], w[2], w[3], w[4] = rk_4ordem(h, 4, t, y0, f)

        for i in 4:N-1
            w[i+1] = w[i] + (h/24) * (55*f(t[i], w[i]) - 59*f(t[i-1], w[i-1]) + 37*f(t[i-2], w[i-2]) - 9*f(t[i-3], w[i-3]))
        end
        
        return w
    end

    # Adams-Moulton 4 ordem 3 passos
    function adams_moulton_4(h, N, t, y0, f::Function)
        w = adams_bashforth_4(h, N, t, y0, f) # predictor

        for i in 3:N-1
            w[i+1] = w[i] + (h/24) * (9*f(t[i+1], w[i+1]) + 19*f(t[i], w[i]) - 5*f(t[i-1], w[i-1]) + f(t[i-2], w[i-2]))
        end
        
        return w
    end

    # SISTEMA DE EQUACOES
    function euler_sistema(h, N, t, y1_0, y2_0, f1::Function, f2::Function)
        w1 = zeros(Float64, N)
        w1[1] = y1_0

        w2 = zeros(Float64, N)
        w2[1] = y2_0

        for i=1:N-1
            w1[i+1] = w1[i] + h*f1(t[i], w1[i], w2[i])
            w2[i+1] = w2[i] + h*f2(t[i], w1[i], w2[i])
        end

        return w1,w2
    end

    function rk_4ordem_sistema(h, N, t, y1_0, y2_0, f1::Function, f2::Function)
        w1 = zeros(Float64, N)
        w1[1] = y1_0

        w2 = zeros(Float64, N)
        w2[1] = y2_0

        for i=1:N-1
            k1 = f1(t[i], w1[i], w2[i])
            k2 = f1(t[i]+h/2, w1[i]+k1*h/2, w2[i]+k1*h/2)
            k3 = f1(t[i]+h/2, w1[i]+k2*h/2, w2[i]+k2*h/2)
            k4 = f1(t[i]+h, w1[i]+k3*h, w2[i]+k2*h/2)
            w1[i+1] = w1[i] + h*(k1+2*k2+2*k3+k4)/6

            k1 = f2(t[i], w1[i], w2[i])
            k2 = f2(t[i]+h/2, w1[i]+k1*h/2, w2[i]+k1*h/2)
            k3 = f2(t[i]+h/2, w1[i]+k2*h/2, w2[i]+k2*h/2)
            k4 = f2(t[i]+h, w1[i]+k3*h, w2[i]+k3*h)
            w2[i+1] = w2[i] + h*(k1+2*k2+2*k3+k4)/6
        end

        return w1,w2
    end

    # Adams-Bashforth 4 ordem 4 passos
    function adams_bashforth_4_sistema(h, N, t, y1_0, y2_0, f1::Function, f2::Function)
        w1 = zeros(Float64, N)
        w2 = zeros(Float64, N)
        aux = rk_4ordem_sistema(h, 4, t, y1_0, y2_0, f1, f2)

        for i in 1:4
            w1[i] = aux[1][i]
            w2[i] = aux[2][i]
        end

        for i in 4:N-1
            w1[i+1] = w1[i] + (h/24) * (55*f1(t[i], w1[i], w2[i]) - 59*f1(t[i-1], w1[i-1], w2[i-1]) + 37*f1(t[i-2], w1[i-2], w2[i-2]) - 9*f1(t[i-3], w1[i-3], w2[i-3]))
            w2[i+1] = w2[i] + (h/24) * (55*f2(t[i], w1[i], w2[i]) - 59*f2(t[i-1], w1[i-1], w2[i-1]) + 37*f2(t[i-2], w1[i-2], w2[i-2]) - 9*f2(t[i-3], w1[i-3], w2[i-3]))
        end
        
        return w1,w2
    end

    # Adams-Moulton 4 ordem 3 passos
    function adams_moulton_4_sistema(h, N, t, y1_0, y2_0, f1::Function, f2::Function)
        w1,w2 = adams_bashforth_4_sistema(h, N, t, y1_0, y2_0, f1, f2) # predictor

        for i in 3:N-1
            w1[i+1] = w1[i] + (h/24) * (9*f1(t[i+1], w1[i+1], w2[i+1]) + 19*f1(t[i], w1[i], w2[i]) - 5*f1(t[i-1], w1[i-1], w2[i-1]) + f1(t[i-2], w1[i-2], w2[i-2]))
            w2[i+1] = w2[i] + (h/24) * (9*f2(t[i+1], w1[i+1], w2[i+1]) + 19*f2(t[i], w1[i], w2[i]) - 5*f2(t[i-1], w1[i-1], w2[i-1]) + f2(t[i-2], w1[i-2], w2[i-2]))
        end
        
        return w1,w2
    end

end
