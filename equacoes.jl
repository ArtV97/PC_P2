# Equacao do Burden
module Burden
    export name, ta, tb, y0, f, fy

    name = "Burden"
    ta = 0
    tb = 2
    y0 = 0.5
    
    # Derivada
    function f(t,y)
        y = y - (t*t) + 1
        return y
    end

    # Solucao analitica
    function fy(t)
        y = (t+1)*(t+1) -0.5*exp(t)
        return y
    end
end

# Equacao da Wikipedia
module Wikipedia
    export name, ta, tb, y0, f, fy

    name = "Wikipedia"
    ta = 0
    tb = 4
    y0 = 1

    # Derivada
    function f(t,y)
        return y
    end

    # Solucao analitica
    function fy(t)
        y = MathConstants.e^t
        return y
    end
end

# Equacao do Chapra
module Chapra
    export name, ta, tb, y0, f, fy

    name = "Chapra"
    ta = 0
    tb = 4
    y0 = 1

    # Derivada
    function f(t,y)
        y = -2 * t^3 + 12 * t^2 - 20 * t + 8.5
        return y
    end

    # Solucao analitica
    function fy(t)
        y = -0.5 * t^4 + 4 * t^3 - 10 * t^2 + 8.5 * t + 1
        return y
    end
end

module Questao3
    export name, y1_0, y2_0, f1, f2, fy

    name = "Questao3"
    ta = 0
    y1_0 = 0
    y2_0 = 0

    # sistema de equacoes
    # y1' =  y2
    # y2' = 1 - y1
    function f1(t, y1, y2)
        return y2
    end
    
    function f2(t, y1, y2)
        return 1 - y1
    end

    # solucao analitica
    function fy(t)
        return 1 - Base.cos(t)
    end
end