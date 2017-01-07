% !! not finished
%position = relative longitudinal position
% e.g. 0.5 is the middle of the segment

function plot_biaxial_test(position, name)

load(name);

z_position = z_0*position
% find the node closest to the position
min_val = 1000;
z_node =1;

for i=1:n_node
    val= (Data(i,1)-z_position)^2;
    if val < min_val
        min_val = val;
        z_node = i;
    end
end

if z_node ==1
    z_node =2;
elseif z_node ==n_node
    z_node = n_node-1;
end

z_TD_L1=0.5*(TD_L1(floor((z_node-1)/(npe-1))*n_Gauss_pt, :) + TD_L1(floor((z_node-1)/(npe-1))*n_Gauss_pt+1, :));
z_TD_L2=0.5*(TD_L2(floor((z_node-1)/(npe-1))*n_Gauss_pt, :) + TD_L2(floor((z_node-1)/(npe-1))*n_Gauss_pt+1, :));
z_TD_mc1=0.5*(TD_mc1(floor((z_node-1)/(npe-1))*n_Gauss_pt, :) + TD_mc1(floor((z_node-1)/(npe-1))*n_Gauss_pt+1, :));
z_TD_mc2=0.5*(TD_mc2(floor((z_node-1)/(npe-1))*n_Gauss_pt, :) + TD_mc2(floor((z_node-1)/(npe-1))*n_Gauss_pt+1, :));
z_TD_mc3=0.5*(TD_mc3(floor((z_node-1)/(npe-1))*n_Gauss_pt, :) + TD_mc3(floor((z_node-1)/(npe-1))*n_Gauss_pt+1, :));
z_TD_mm =0.5*(TD_mm(floor((z_node-1)/(npe-1))*n_Gauss_pt, :) + TD_mm(floor((z_node-1)/(npe-1))*n_Gauss_pt+1, :));

z_TD_L2
itnum=0;
itr_error = 10.0;  %initial value of iteration error
P_i = [0.92; 0.92]; 


while (  itr_error > tol/100)
    itnum = itnum+1;

    %initializing the matrix

    K_ij = [0, 0;0, 0];
    t_ij = [0; 0];


    L1_s=z_TD_L1(n_time+1);                  %lambda_1 (circumflex)
    L2_s=z_TD_L2(n_time+1);
    
    L1 = L1_s * P_i(1);
    L2 = L2_s * P_i(2);

    if n_time+1 <= num_DL  %If there still remains initial mass
        % phi = phi0
        L_k3=sqrt( (L1*sin(phi0))^2+(L2*cos(phi0))^2);
        L_k4=sqrt( (L1*sin(pi-phi0))^2+(L2*cos(pi-phi0))^2);

        dL_k3dL1=L1*sin(phi0)^2 / L_k3;  %dL_k4dL1=dL_k3dL1
        dL_k3dL2=L2*cos(phi0)^2 / L_k3;
        ddL_k3ddL1=(L2*sin(phi0)*cos(phi0))^2 /L_k3^3;
        ddL_k3ddL2=(L1*sin(phi0)*cos(phi0))^2 /L_k3^3;
        ddL_k3dL1dL2= - L1*L2*(sin(phi0)*cos(phi0))^2 / L_k3^3;
        dL_k4dL1=L1*sin(pi-phi0)^2 / L_k4;
        dL_k4dL2=L2*cos(pi-phi0)^2 / L_k4;
        ddL_k4ddL1=(L2*sin(pi-phi0)*cos(pi-phi0))^2 /L_k4^3;
        ddL_k4ddL2=(L1*sin(pi-phi0)*cos(pi-phi0))^2 /L_k4^3;
        ddL_k4dL1dL2= - L1*L2*(sin(pi-phi0)*cos(pi-phi0))^2 / L_k4^3;

        Ln_k1 = G_h*L1;
        Ln_k2 = G_h*L2;
        Ln_k3 = G_h*L_k3;
        Ln_k4 = G_h*L_k4;

        Mf_r=Mf0*DQ_c(n_time+1);
        Mt = Mf_r(1)+Mf_r(2)+Mf_r(3)+Mf_r(4);

        dPsi_k1dL1 = G_h*dWcdx(Ln_k1);
        dPsi_k3dL1 = G_h*dWcdx(Ln_k3)*dL_k3dL1;
        dPsi_k4dL1 = G_h*dWcdx(Ln_k4)*dL_k4dL1;
        dPsi_k2dL2 = G_h*dWcdx(Ln_k2);
        dPsi_k3dL2 = G_h*dWcdx(Ln_k3)*dL_k3dL2;
        dPsi_k4dL2 = G_h*dWcdx(Ln_k4)*dL_k4dL2;
        
        ddPsi_k1ddL1 = G_h^2*ddWcddx(Ln_k1);
        ddPsi_k3ddL1 = G_h^2*ddWcddx(Ln_k3)*dL_k3dL1^2+G_h*dWcdx(Ln_k3)*ddL_k3ddL1;
        ddPsi_k4ddL1 = G_h^2*ddWcddx(Ln_k4)*dL_k4dL1^2+G_h*dWcdx(Ln_k4)*ddL_k4ddL1;
        ddPsi_k2ddL2 = G_h^2*ddWcddx(Ln_k2);
        ddPsi_k3ddL2 = G_h^2*ddWcddx(Ln_k3)*dL_k3dL2^2+G_h*dWcdx(Ln_k3)*ddL_k3ddL2;
        ddPsi_k4ddL2 = G_h^2*ddWcddx(Ln_k4)*dL_k4dL2^2+G_h*dWcdx(Ln_k4)*ddL_k4ddL2;
        ddPsi_k3dL1dL2 = G_h^2*ddWcddx(Ln_k3)*dL_k3dL2*dL_k3dL1+G_h*dWcdx(Ln_k3)*ddL_k3dL1dL2;
        ddPsi_k4dL1dL2 = G_h^2*ddWcddx(Ln_k4)*dL_k4dL2*dL_k4dL1+G_h*dWcdx(Ln_k4)*ddL_k4dL1dL2;;
        
        Mm = Mm0*DQ_m(n_time+1)*Rm_exp(z_position, z_0, n_time*dt, init_dmg_t);
        Mt = Mt+Mm;
        Ln_m = G_m * L1;
        dPsi_mdL1 = G_m*dWmdx(Ln_m);
        ddPsi_mddL1 = G_m^2*ddWmddx(Ln_m);
        
        K_ij(1,1) = K_ij(1,1)+  Mf_r(1)*L1_s*(dPsi_k1dL1+L1*ddPsi_k1ddL1) ...
            + Mf_r(3)*L1_s*(dPsi_k3dL1+L1*ddPsi_k3ddL1)...
            + Mf_r(4)*L1_s*(dPsi_k4dL1+L1*ddPsi_k4ddL1)...
            + Mm*L1_s*(dPsi_mdL1+L1*ddPsi_k1ddL1);
        K_ij(1,2) = K_ij(1,2)+ Mf_r(3)*L2_s*L1*ddPsi_k3dL1dL2...
            + Mf_r(4)*L2_s*L1*ddPsi_k4dL1dL2;
        K_ij(2,1) = K_ij(2,1)+ Mf_r(3)*L1_s*L2*ddPsi_k3dL1dL2...
            + Mf_r(4)*L1_s*L2*ddPsi_k4dL1dL2;
        K_ij(2,2) = K_ij(2,2)+  Mf_r(2)*L2_s*(dPsi_k2dL2+L2*ddPsi_k2ddL2) ...
            + Mf_r(3)*L2_s*(dPsi_k3dL2+L2*ddPsi_k3ddL2)...
            + Mf_r(4)*L2_s*(dPsi_k4dL2+L2*ddPsi_k4ddL2);
        
        t_ij(1) = t_ij(1)+ Mf_r(1)*L1*dPsi_k1dL1 ...
            + Mf_r(3)*L1*dPsi_k3dL1...
            + Mf_r(4)*L1*dPsi_k4dL1...
            + Mm*L1*dPsi_mdL1;
        t_ij(2) = t_ij(2)+ Mf_r(2)*L2*dPsi_k2dL2 ...
            + Mf_r(3)*L2*dPsi_k3dL2...
            + Mf_r(4)*L2*dPsi_k4dL2;
        % we assume that the initial SM removed with the elastin
        % degradation
    else
        Mt =0.0;
        Mm =0.0;
    end

    % strain energy due to elastin layer

    Me = Me0*Rm_exp(z_position, z_0, n_time*dt, init_dmg_t);
    Mt = Mt+Me;

    Ln_e1 = G_e(1)*L1;
    Ln_e2 = G_e(2)*L2;

    dPsi_edL1= G_e(1)*dWedx(1,Ln_e1, Ln_e2);
    dPsi_edL2= G_e(2)*dWedx(2,Ln_e1, Ln_e2);
    ddPsi_eddL1=(G_e(1))^2*ddWeddx(1,Ln_e1,Ln_e2);
    ddPsi_eddL2=(G_e(2))^2*ddWeddx(2,Ln_e1,Ln_e2);
    ddPsi_edL1dL2 = G_e(1)*G_e(2)*ddWeddx(3,Ln_e1,Ln_e2);

    t_ij(1) = t_ij(1)+Me*L1*dPsi_edL1;
    t_ij(2) = t_ij(2)+Me*L2*dPsi_edL2;
    
    K_ij(1,1) = K_ij(1,1)+Me*L1_s*(dPsi_edL1+L1*ddPsi_eddL1);
    K_ij(1,2)= K_ij(1,2)+Me*L2_s*L1*ddPsi_edL1dL2;
    K_ij(2,1)= K_ij(2,1)+Me*L1_s*L2*ddPsi_edL1dL2;
    K_ij(2,2) = K_ij(2,2)+Me*L2_s*(dPsi_edL2+L2*ddPsi_eddL2);

    %-----------------------------------------------------------
    %       numerical integration
    % int(0 to t) [ ff(tau)]dtau = ff(0)*0.5*dtau+
    % ff(dtau)*dtau+ff(2*dtau)*dtau+ . .+ff(t)*0.5*dtau
    %-----------------------------------------------------------
    if n_time+1 <= num_DL 
        n_tau0 = 1;
    else
        n_tau0= n_time+1-num_DL;
    end
    for n_tau=n_tau0: n_time+1
        tau=(n_tau-1)*dt;
        L1_tau=z_TD_L1(n_tau);
        L2_tau=z_TD_L2(n_tau);

        phi = phi0;  % phi(\tau) at the reference configuration

        L_k3=sqrt( (L1*sin(phi))^2+(L2*cos(phi))^2);
        L_k4=sqrt( (L1*sin(pi-phi))^2+(L2*cos(pi-phi))^2);

        dL_k3dL1=L1*sin(phi)^2 / L_k3;  %dL_k4dL1=dL_k3dL1
        dL_k3dL2=L2*cos(phi)^2 / L_k3;
        ddL_k3ddL1=(L2*sin(phi)*cos(phi))^2 /L_k3^3;
        ddL_k3ddL2=(L1*sin(phi)*cos(phi))^2 /L_k3^3;
        ddL_k3dL1dL2= - L1*L2*(sin(phi)*cos(phi))^2 / L_k3^3;
        dL_k4dL1=L1*sin(pi-phi)^2 / L_k4;
        dL_k4dL2=L2*cos(pi-phi)^2 / L_k4;
        ddL_k4ddL1=(L2*sin(pi-phi)*cos(pi-phi))^2 /L_k4^3;
        ddL_k4ddL2=(L1*sin(pi-phi)*cos(pi-phi))^2 /L_k4^3;
        ddL_k4dL1dL2= - L1*L2*(sin(pi-phi)*cos(pi-phi))^2 / L_k4^3;

        L_k3_tau=sqrt( (L1_tau*sin(phi))^2+(L2_tau*cos(phi))^2);
        L_k4_tau = L_k3_tau;

        if ((n_tau==n_tau0) || (n_tau==nt+1))
            w_dt=dt*0.5;
        else
            w_dt=dt;
        end

        mc_tau=[TD_mc1(count,n_tau), TD_mc2(count, n_tau), ...
            TD_mc3(count, n_tau), TD_mc3(count, n_tau)]*q_i(0, n_time*dt-tau);
        mm_tau = TD_mm(count, n_tau) * q_i(1, n_time*dt-tau)...
            *Rm_exp(z_position, z_0, n_time*dt, init_dmg_t);
        Mm = Mm+mm_tau*w_dt;
        Mt = Mt+ (mc_tau(1)+mc_tau(2)+mc_tau(3)+mc_tau(4))*w_dt...
            +mm_tau*w_dt;

        Ln_k1 = G_h*L1 / L1_tau;
        Ln_k2 = G_h*L2 / L2_tau;
        Ln_k3 = G_h*L_k3 / L_k3_tau;
        Ln_k4 = G_h*L_k4 / L_k4_tau;
        Ln_m =  G_m*L1/L1_tau;

        dPsi_k1dL1 = (G_h/L1_tau)*dWcdx(Ln_k1);
        dPsi_k3dL1 = (G_h/L_k3_tau)*dWcdx(Ln_k3)*dL_k3dL1;
        dPsi_k4dL1 = (G_h/L_k4_tau)*dWcdx(Ln_k4)*dL_k4dL1;
        dPsi_k2dL2 = (G_h/L2_tau)*dWcdx(Ln_k2);
        dPsi_k3dL2 = (G_h/L_k3_tau)*dWcdx(Ln_k3)*dL_k3dL2;
        dPsi_k4dL2 = (G_h/L_k4_tau)*dWcdx(Ln_k4)*dL_k4dL2;

        ddPsi_k1ddL1 = (G_h/L1_tau)^2*ddWcddx(Ln_k1);
        ddPsi_k3ddL1 = (G_h/L_k3_tau)^2*ddWcddx(Ln_k3)*dL_k3dL1^2+(G_h/L_k3_tau)*dWcdx(Ln_k3)*ddL_k3ddL1;
        ddPsi_k4ddL1 = (G_h/L_k4_tau)^2*ddWcddx(Ln_k4)*dL_k4dL1^2+(G_h/L_k4_tau)*dWcdx(Ln_k4)*ddL_k4ddL1;
        ddPsi_k2ddL2 = (G_h/L2_tau)^2*ddWcddx(Ln_k2);
        ddPsi_k3ddL2 = (G_h/L_k3_tau)^2*ddWcddx(Ln_k3)*dL_k3dL2^2+(G_h/L_k3_tau)*dWcdx(Ln_k3)*ddL_k3ddL2;
        ddPsi_k4ddL2 = (G_h/L_k4_tau)^2*ddWcddx(Ln_k4)*dL_k4dL2^2+(G_h/L_k4_tau)*dWcdx(Ln_k4)*ddL_k4ddL2;
        ddPsi_k3dL1dL2 = (G_h/L_k3_tau)^2*ddWcddx(Ln_k3)*dL_k3dL2*dL_k3dL1+(G_h/L_k3_tau)*dWcdx(Ln_k3)*ddL_k3dL1dL2;
        ddPsi_k4dL1dL2 = (G_h/L_k4_tau)^2*ddWcddx(Ln_k4)*dL_k4dL2*dL_k4dL1+(G_h/L_k4_tau)*dWcdx(Ln_k4)*ddL_k4dL1dL2;

        dPsi_mdL1 = (G_m/L1_tau)*dWmdx(Ln_m);
        ddPsi_mddL1 = (G_m/L1_tau)^2*ddWmddx(Ln_m);

        K_ij(1,1) = K_ij(1,1)+  (mc_tau(1)*L1_s*(dPsi_k1dL1+L1*ddPsi_k1ddL1) ...
            + mc_tau(3)*L1_s*(dPsi_k3dL1+L1*ddPsi_k3ddL1)...
            + mc_tau(4)*L1_s*(dPsi_k4dL1+L1*ddPsi_k4ddL1)...
            + mm_tau*L1_s*(dPsi_mdL1+L1*ddPsi_k1ddL1))*w_dt;
        K_ij(1,2) = K_ij(1,2)+ (mc_tau(3)*L2_s*L1*ddPsi_k3dL1dL2...
            + mc_tau(4)*L2_s*L1*ddPsi_k4dL1dL2)*w_dt;
        K_ij(2,1) = K_ij(2,1)+ (mc_tau(3)*L1_s*L2*ddPsi_k3dL1dL2...
            + mc_tau(4)*L1_s*L2*ddPsi_k4dL1dL2)*w_dt;
        K_ij(2,2) = K_ij(2,2)+  (mc_tau(2)*L2_s*(dPsi_k2dL2+L2*ddPsi_k2ddL2) ...
            + mc_tau(3)*L2_s*(dPsi_k3dL2+L2*ddPsi_k3ddL2)...
            + mc_tau(4)*L2_s*(dPsi_k4dL2+L2*ddPsi_k4ddL2))*w_dt;

        t_ij(1) = t_ij(1)+ (mc_tau(1)*L1*dPsi_k1dL1 ...
            + mc_tau(3)*L1*dPsi_k3dL1...
            + mc_tau(4)*L1*dPsi_k4dL1...
            + mm_tau*L1*dPsi_mdL1)*w_dt;
        t_ij(2) = t_ij(2)+ (mc_tau(2)*L2*dPsi_k2dL2 ...
            + mc_tau(3)*L2*dPsi_k3dL2...
            + mc_tau(4)*L2*dPsi_k4dL2)*w_dt;
    end


    %glf
    %glk
    %solving matrix

    dP_i=K_ij\t_ij;

    P_i=P_i-dP_i;
    % error
    itr_error = sum(dP_i.^2)/sum(P_i.^2);
end
 

P_i0=P_i

t_11=0.0;
t_22=0.0;
n=0;
while (t_11<=7e+4 & t_22<=7e+4)
    t_11 = n*700;
    t_22 = n*700;
    n=n+1;
    itnum=0;
    itr_error = 10.0;  %initial value of iteration error
    P_i =P_i0; %initial guess
    while ( itr_error > tol/100)
        itnum = itnum+1;

        %initializing the matrix

        K_ij = [0, 0;0, 0];
        t_ij = [0; 0];


        L1_s=z_TD_L1(n_time+1);                  %lambda_1 (circumflex)
        L2_s=z_TD_L2(n_time+1);

        L1 = L1_s * P_i(1);
        L2 = L2_s * P_i(2);

        if n_time+1 <= num_DL  %If there still remains initial mass
            % phi = phi0
            L_k3=sqrt( (L1*sin(phi0))^2+(L2*cos(phi0))^2);
            L_k4=sqrt( (L1*sin(pi-phi0))^2+(L2*cos(pi-phi0))^2);

            dL_k3dL1=L1*sin(phi0)^2 / L_k3;  %dL_k4dL1=dL_k3dL1
            dL_k3dL2=L2*cos(phi0)^2 / L_k3;
            ddL_k3ddL1=(L2*sin(phi0)*cos(phi0))^2 /L_k3^3;
            ddL_k3ddL2=(L1*sin(phi0)*cos(phi0))^2 /L_k3^3;
            ddL_k3dL1dL2= - L1*L2*(sin(phi0)*cos(phi0))^2 / L_k3^3;
            dL_k4dL1=L1*sin(pi-phi0)^2 / L_k4;
            dL_k4dL2=L2*cos(pi-phi0)^2 / L_k4;
            ddL_k4ddL1=(L2*sin(pi-phi0)*cos(pi-phi0))^2 /L_k4^3;
            ddL_k4ddL2=(L1*sin(pi-phi0)*cos(pi-phi0))^2 /L_k4^3;
            ddL_k4dL1dL2= - L1*L2*(sin(pi-phi0)*cos(pi-phi0))^2 / L_k4^3;

            Ln_k1 = G_h*L1;
            Ln_k2 = G_h*L2;
            Ln_k3 = G_h*L_k3;
            Ln_k4 = G_h*L_k4;

            Mf_r=Mf0*DQ_c(n_time+1);
            Mt = Mf_r(1)+Mf_r(2)+Mf_r(3)+Mf_r(4);

            dPsi_k1dL1 = G_h*dWcdx(Ln_k1);
            dPsi_k3dL1 = G_h*dWcdx(Ln_k3)*dL_k3dL1;
            dPsi_k4dL1 = G_h*dWcdx(Ln_k4)*dL_k4dL1;
            dPsi_k2dL2 = G_h*dWcdx(Ln_k2);
            dPsi_k3dL2 = G_h*dWcdx(Ln_k3)*dL_k3dL2;
            dPsi_k4dL2 = G_h*dWcdx(Ln_k4)*dL_k4dL2;

            ddPsi_k1ddL1 = G_h^2*ddWcddx(Ln_k1);
            ddPsi_k3ddL1 = G_h^2*ddWcddx(Ln_k3)*dL_k3dL1^2+G_h*dWcdx(Ln_k3)*ddL_k3ddL1;
            ddPsi_k4ddL1 = G_h^2*ddWcddx(Ln_k4)*dL_k4dL1^2+G_h*dWcdx(Ln_k4)*ddL_k4ddL1;
            ddPsi_k2ddL2 = G_h^2*ddWcddx(Ln_k2);
            ddPsi_k3ddL2 = G_h^2*ddWcddx(Ln_k3)*dL_k3dL2^2+G_h*dWcdx(Ln_k3)*ddL_k3ddL2;
            ddPsi_k4ddL2 = G_h^2*ddWcddx(Ln_k4)*dL_k4dL2^2+G_h*dWcdx(Ln_k4)*ddL_k4ddL2;
            ddPsi_k3dL1dL2 = G_h^2*ddWcddx(Ln_k3)*dL_k3dL2*dL_k3dL1+G_h*dWcdx(Ln_k3)*ddL_k3dL1dL2;
            ddPsi_k4dL1dL2 = G_h^2*ddWcddx(Ln_k4)*dL_k4dL2*dL_k4dL1+G_h*dWcdx(Ln_k4)*ddL_k4dL1dL2;;

            Mm = Mm0*DQ_m(n_time+1)*Rm_exp(z_position, z_0, n_time*dt, init_dmg_t);
            Mt = Mt+Mm;
            Ln_m = G_m * L1;
            dPsi_mdL1 = G_m*dWmdx(Ln_m);
            ddPsi_mddL1 = G_m^2*ddWmddx(Ln_m);

            K_ij(1,1) = K_ij(1,1)+  Mf_r(1)*L1_s*(dPsi_k1dL1+L1*ddPsi_k1ddL1) ...
                + Mf_r(3)*L1_s*(dPsi_k3dL1+L1*ddPsi_k3ddL1)...
                + Mf_r(4)*L1_s*(dPsi_k4dL1+L1*ddPsi_k4ddL1)...
                + Mm*L1_s*(dPsi_mdL1+L1*ddPsi_k1ddL1);
            K_ij(1,2) = K_ij(1,2)+ Mf_r(3)*L2_s*L1*ddPsi_k3dL1dL2...
                + Mf_r(4)*L2_s*L1*ddPsi_k4dL1dL2;
            K_ij(2,1) = K_ij(2,1)+ Mf_r(3)*L1_s*L2*ddPsi_k3dL1dL2...
                + Mf_r(4)*L1_s*L2*ddPsi_k4dL1dL2;
            K_ij(2,2) = K_ij(2,2)+  Mf_r(2)*L2_s*(dPsi_k2dL2+L2*ddPsi_k2ddL2) ...
                + Mf_r(3)*L2_s*(dPsi_k3dL2+L2*ddPsi_k3ddL2)...
                + Mf_r(4)*L2_s*(dPsi_k4dL2+L2*ddPsi_k4ddL2);

            t_ij(1) = t_ij(1)+ Mf_r(1)*L1*dPsi_k1dL1 ...
                + Mf_r(3)*L1*dPsi_k3dL1...
                + Mf_r(4)*L1*dPsi_k4dL1...
                + Mm*L1*dPsi_mdL1;
            t_ij(2) = t_ij(2)+ Mf_r(2)*L2*dPsi_k2dL2 ...
                + Mf_r(3)*L2*dPsi_k3dL2...
                + Mf_r(4)*L2*dPsi_k4dL2;
            % we assume that the initial SM removed with the elastin
            % degradation
        else
            Mt =0.0;
            Mm =0.0;
        end

        % strain energy due to elastin layer

        Me = Me0*Rm_exp(z_position, z_0, n_time*dt, init_dmg_t);
        Mt = Mt+Me;

        Ln_e1 = G_e(1)*L1;
        Ln_e2 = G_e(2)*L2;

        dPsi_edL1= G_e(1)*dWedx(1,Ln_e1, Ln_e2);
        dPsi_edL2= G_e(2)*dWedx(2,Ln_e1, Ln_e2);
        ddPsi_eddL1=(G_e(1))^2*ddWeddx(1,Ln_e1,Ln_e2);
        ddPsi_eddL2=(G_e(2))^2*ddWeddx(2,Ln_e1,Ln_e2);
        ddPsi_edL1dL2 = G_e(1)*G_e(2)*ddWeddx(3,Ln_e1,Ln_e2);

        t_ij(1) = t_ij(1)+Me*L1*dPsi_edL1;
        t_ij(2) = t_ij(2)+Me*L2*dPsi_edL2;

        K_ij(1,1) = K_ij(1,1)+Me*L1_s*(dPsi_edL1+L1*ddPsi_eddL1);
        K_ij(1,2)= K_ij(1,2)+Me*L2_s*L1*ddPsi_edL1dL2;
        K_ij(2,1)= K_ij(2,1)+Me*L1_s*L2*ddPsi_edL1dL2;
        K_ij(2,2) = K_ij(2,2)+Me*L2_s*(dPsi_edL2+L2*ddPsi_eddL2);

        %-----------------------------------------------------------
        %       numerical integration
        % int(0 to t) [ ff(tau)]dtau = ff(0)*0.5*dtau+
        % ff(dtau)*dtau+ff(2*dtau)*dtau+ . .+ff(t)*0.5*dtau
        %-----------------------------------------------------------
        if n_time+1 <= num_DL
            n_tau0 = 1;
        else
            n_tau0= n_time+1-num_DL;
        end
        for n_tau=n_tau0: n_time+1
            tau=(n_tau-1)*dt;
            L1_tau=z_TD_L1(n_tau);
            L2_tau=z_TD_L2(n_tau);

            phi = phi0;  % phi(\tau) at the reference configuration

            L_k3=sqrt( (L1*sin(phi))^2+(L2*cos(phi))^2);
            L_k4=sqrt( (L1*sin(pi-phi))^2+(L2*cos(pi-phi))^2);

            dL_k3dL1=L1*sin(phi)^2 / L_k3;  %dL_k4dL1=dL_k3dL1
            dL_k3dL2=L2*cos(phi)^2 / L_k3;
            ddL_k3ddL1=(L2*sin(phi)*cos(phi))^2 /L_k3^3;
            ddL_k3ddL2=(L1*sin(phi)*cos(phi))^2 /L_k3^3;
            ddL_k3dL1dL2= - L1*L2*(sin(phi)*cos(phi))^2 / L_k3^3;
            dL_k4dL1=L1*sin(pi-phi)^2 / L_k4;
            dL_k4dL2=L2*cos(pi-phi)^2 / L_k4;
            ddL_k4ddL1=(L2*sin(pi-phi)*cos(pi-phi))^2 /L_k4^3;
            ddL_k4ddL2=(L1*sin(pi-phi)*cos(pi-phi))^2 /L_k4^3;
            ddL_k4dL1dL2= - L1*L2*(sin(pi-phi)*cos(pi-phi))^2 / L_k4^3;

            L_k3_tau=sqrt( (L1_tau*sin(phi))^2+(L2_tau*cos(phi))^2);
            L_k4_tau = L_k3_tau;

            if ((n_tau==n_tau0) || (n_tau==nt+1))
                w_dt=dt*0.5;
            else
                w_dt=dt;
            end

            mc_tau=[TD_mc1(count,n_tau), TD_mc2(count, n_tau), ...
                TD_mc3(count, n_tau), TD_mc3(count, n_tau)]*q_i(0, n_time*dt-tau);
            mm_tau = TD_mm(count, n_tau) * q_i(1, n_time*dt-tau)...
                *Rm_exp(z_position, z_0, n_time*dt, init_dmg_t);
            Mm = Mm+mm_tau*w_dt;
            Mt = Mt+ (mc_tau(1)+mc_tau(2)+mc_tau(3)+mc_tau(4))*w_dt...
                +mm_tau*w_dt;

            Ln_k1 = G_h*L1 / L1_tau;
            Ln_k2 = G_h*L2 / L2_tau;
            Ln_k3 = G_h*L_k3 / L_k3_tau;
            Ln_k4 = G_h*L_k4 / L_k4_tau;
            Ln_m =  G_m*L1/L1_tau;

            dPsi_k1dL1 = (G_h/L1_tau)*dWcdx(Ln_k1);
            dPsi_k3dL1 = (G_h/L_k3_tau)*dWcdx(Ln_k3)*dL_k3dL1;
            dPsi_k4dL1 = (G_h/L_k4_tau)*dWcdx(Ln_k4)*dL_k4dL1;
            dPsi_k2dL2 = (G_h/L2_tau)*dWcdx(Ln_k2);
            dPsi_k3dL2 = (G_h/L_k3_tau)*dWcdx(Ln_k3)*dL_k3dL2;
            dPsi_k4dL2 = (G_h/L_k4_tau)*dWcdx(Ln_k4)*dL_k4dL2;

            ddPsi_k1ddL1 = (G_h/L1_tau)^2*ddWcddx(Ln_k1);
            ddPsi_k3ddL1 = (G_h/L_k3_tau)^2*ddWcddx(Ln_k3)*dL_k3dL1^2+(G_h/L_k3_tau)*dWcdx(Ln_k3)*ddL_k3ddL1;
            ddPsi_k4ddL1 = (G_h/L_k4_tau)^2*ddWcddx(Ln_k4)*dL_k4dL1^2+(G_h/L_k4_tau)*dWcdx(Ln_k4)*ddL_k4ddL1;
            ddPsi_k2ddL2 = (G_h/L2_tau)^2*ddWcddx(Ln_k2);
            ddPsi_k3ddL2 = (G_h/L_k3_tau)^2*ddWcddx(Ln_k3)*dL_k3dL2^2+(G_h/L_k3_tau)*dWcdx(Ln_k3)*ddL_k3ddL2;
            ddPsi_k4ddL2 = (G_h/L_k4_tau)^2*ddWcddx(Ln_k4)*dL_k4dL2^2+(G_h/L_k4_tau)*dWcdx(Ln_k4)*ddL_k4ddL2;
            ddPsi_k3dL1dL2 = (G_h/L_k3_tau)^2*ddWcddx(Ln_k3)*dL_k3dL2*dL_k3dL1+(G_h/L_k3_tau)*dWcdx(Ln_k3)*ddL_k3dL1dL2;
            ddPsi_k4dL1dL2 = (G_h/L_k4_tau)^2*ddWcddx(Ln_k4)*dL_k4dL2*dL_k4dL1+(G_h/L_k4_tau)*dWcdx(Ln_k4)*ddL_k4dL1dL2;

            dPsi_mdL1 = (G_m/L1_tau)*dWmdx(Ln_m);
            ddPsi_mddL1 = (G_m/L1_tau)^2*ddWmddx(Ln_m);

            K_ij(1,1) = K_ij(1,1)+  (mc_tau(1)*L1_s*(dPsi_k1dL1+L1*ddPsi_k1ddL1) ...
                + mc_tau(3)*L1_s*(dPsi_k3dL1+L1*ddPsi_k3ddL1)...
                + mc_tau(4)*L1_s*(dPsi_k4dL1+L1*ddPsi_k4ddL1)...
                + mm_tau*L1_s*(dPsi_mdL1+L1*ddPsi_k1ddL1))*w_dt;
            K_ij(1,2) = K_ij(1,2)+ (mc_tau(3)*L2_s*L1*ddPsi_k3dL1dL2...
                + mc_tau(4)*L2_s*L1*ddPsi_k4dL1dL2)*w_dt;
            K_ij(2,1) = K_ij(2,1)+ (mc_tau(3)*L1_s*L2*ddPsi_k3dL1dL2...
                + mc_tau(4)*L1_s*L2*ddPsi_k4dL1dL2)*w_dt;
            K_ij(2,2) = K_ij(2,2)+  (mc_tau(2)*L2_s*(dPsi_k2dL2+L2*ddPsi_k2ddL2) ...
                + mc_tau(3)*L2_s*(dPsi_k3dL2+L2*ddPsi_k3ddL2)...
                + mc_tau(4)*L2_s*(dPsi_k4dL2+L2*ddPsi_k4ddL2))*w_dt;

            t_ij(1) = t_ij(1)+ (mc_tau(1)*L1*dPsi_k1dL1 ...
                + mc_tau(3)*L1*dPsi_k3dL1...
                + mc_tau(4)*L1*dPsi_k4dL1...
                + mm_tau*L1*dPsi_mdL1)*w_dt;
            t_ij(2) = t_ij(2)+ (mc_tau(2)*L2*dPsi_k2dL2 ...
                + mc_tau(3)*L2*dPsi_k3dL2...
                + mc_tau(4)*L2*dPsi_k4dL2)*w_dt;
        end


        %glf
        %glk
        %solving matrix
        K_ij=rho/Mt*K_ij;
        t_ij=rho/Mt*t_ij-[t_11;t_22];

        dP_i=K_ij\t_ij;

        P_i=P_i-0.5*dP_i;
        % error
        itr_error = sum(dP_i.^2)/sum(P_i.^2);
    end

    plot(P_i(1)/P_i0(1), t_11, 'ro');
    hold on
    plot(P_i(2)/P_i0(2), t_22, 'o');
end
clear *
