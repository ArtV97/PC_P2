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