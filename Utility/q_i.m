%survival function
function y=q_i(cn,t, kq_c, kq_m)

if(cn==0) %collagen: kq_c (1/day),   t (day)
    y= exp(-kq_c*t);
elseif(cn==1)  %SM
    y= exp(-kq_m*t);
elseif(cn==2) %elastin
    y=1;
else
    error('survival function error');
end
