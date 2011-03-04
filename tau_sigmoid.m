function[tau] = tau_sigmoid(a, b, c, d, e, V)
        tau = (a + b/(1+exp(c*(V+d)/e)))*10^-3; %[seconds]
end