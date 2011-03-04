% Sigmoid functions for calculating channel dynamics
function [answer] = sigmoid(a,b,c, V) % a is direction of curve, -b is half (in)activation voltage, V is votlage.
    answer = 1/(1+exp(a*(V+b)/c));
end