% Return the exponential shape of mass

function y=Rm_exp(ID, Z, z_0, t, p_t)

if t <= p_t
    y=1.0;
else
    c1=0.7;
    c2=6.0;
    Zc=z_0/2;
    Z_0=Z-Zc;
    rate= 1/40;
  if ID==0 %elastin & SM
        f0=0.99*exp(-c1*Z_0^c2);
  elseif ID==1 %Collagen
        f0=0.2*exp(-c1*Z_0^c2);
  end
    if f0< 0.01
        f0= 0.01-5*(0.01-f0);
        if f0<0.0
            f0=0.0;
        end
    end
    y=1.0-f0*(1.0-exp(-rate*(t-p_t)));

end
