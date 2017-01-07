function P=P_AAA(t, op_time, P_a)

if t<op_time
    P=P_a;
else
    k=1/40;
    v=exp(-1.0*k*(t-op_time));
    P=P_a*v + P_a*0.4*(1-v);
end

