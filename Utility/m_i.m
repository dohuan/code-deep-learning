%mass production: kinetics of growth

function y=m_i(n_c, m_basal,dn_stress, Kc)

y=n_c*m_basal*(Kc*dn_stress+1.0);
if y<0
    y=0.0;
end