function y=ddWcddx(x, kc)

if x>= 1.0
    dQdx=4.0*kc(3)*x*(x^2-1);

    y=kc(2)*(3*x^2-1+(x^3-x)*dQdx)*exp(kc(3)*(x^2-1)^2);
else
    dQdx=4.0*kc(3)*x*(x^2-1);
 
    y=0.1 * kc(2)*(3*x^2-1+(x^3-x)*dQdx)*exp(kc(3)*(x^2-1)^2);
end

