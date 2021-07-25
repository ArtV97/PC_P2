# Equacao do Burden
module Burden
    export ta, tb, y0, f, fy

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
    export ta, tb, y0, f, fy

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
    export ta, tb, y0, f, fy

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
        y = t^4 + 4 * t^3 - 10 * t^2 + 8.5 * t + 1
        return y
    end
end