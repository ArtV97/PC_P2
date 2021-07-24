# Equacao do Burden

# derivada
function burden_f(t,y)
    y = y - (t*t) + 1
    return y
end

# Solucao analitica Burden
function burden_fy(t)
    y = (t+1)*(t+1) -0.5*exp(t)
    return y
end