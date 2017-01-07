function y=ddWmddx(x, kc)

dQdx=4.0*kc(5)*x*(x^2-1);

y=kc(4)*(3*x^2-1+(x^3-x)*dQdx)*exp(kc(5)*(x^2-1)^2);

