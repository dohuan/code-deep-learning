function y=dWedx(k, L1, L2, kc)

if k==1
   y=kc(1)*L1*(1-1/(L1^4*L2^2));
elseif k==2
   y=kc(1)*L2*(1-1/(L2^4*L1^2));
else
    error('Error in dWedx.m\n');
end

