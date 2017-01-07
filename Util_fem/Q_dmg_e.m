%(JSW) Damage function for elastin (represents loss of elastin due to 
%       stretch beyond the ultimate stretch)
%   Q_dmg_e: survival fraction
%   Ln_e: max current elastin stretch (from natural state)

% NOTE: E=0.5*(L^2-1), where E is strain and L is stretch???

function y=Q_dmg_e(Ln_e1,Ln_e2,Le_low)
%  Q_dmg_e=1-e^[k_dmg*(Ln_e-L_low)]

    I1 = Ln_e1^2 + Ln_e2^2 + 1;
    L_low=Le_low; %low range of elastin ultimate stretch (loss=(1-low_per)%)
    L_high=L_low+1; %high range of elastin ultimate stretch (loss=100%)
    low_frac=0.99; %fraction of elastin remaining at L_low
    
    %Calculate k_dmg, the damage factor
    k_dmg=log(1-low_frac)/(L_low-L_high);
    
    %Calculate Q_dmg_e, the damage survival fraction
    if I1 >= L_high
        y=0; %no survival beyond the high range
    else
        y=1-exp(k_dmg*(I1-L_high));
    end
end