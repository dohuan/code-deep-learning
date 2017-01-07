function y=dWcdx(x, kc)

if x>=1.0
    y=kc(2)*(x*(x^2-1))*exp(kc(3)*(x^2-1)^2);
else
    y=0.1*kc(2)*(x*(x^2-1))*exp(kc(3)*(x^2-1)^2);
end

