function y=ddWeddx(k,L1,L2,kc)

if k==1
   y=kc(1)*( (1-1/(L1^4*L2^2)) + 4/(L1^4*L2^2));
elseif k==2
   y=kc(1)*( (1-1/(L2^4*L1^2)) + L2*4/(L2^5*L1^2));
elseif k==3
    y=kc(1)*2/(L1*L2)^3;
else
    error('Error in dWedx.m\n');
end

