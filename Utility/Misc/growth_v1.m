function growth(days,k_sigma_f, k_sigma_m, name, length, n_dt, print_step)

global kc P_a r_h H_h nu_e0 nu_f0 nu_m0 phi0 G_h G_e G_m ...
    Sa La_M La_0 sigma_f0 sigma_m0 n_elm kq_c kq_m age_max sh type op_time

n_Gauss_pt=3;        % number of gauss points: choose from 1 to 6
type = 'quadratic';  %type of the shape function 'linear' or 'quadratic'
flag_LinearC = 0;    % if flag_LinearC==1, then it calculate the linearized stiffness
bc = 0;  % 0: load-free bounday condition, 1:fixed bc 

r_stent = 1.1;
k_stent = 0.0; %max =20


rho=1050;                         %density of an artery=1050 kg/m^3 
M0=H_h*rho ;   % total mass area density
Me0=M0*nu_e0;
Mf0=M0*nu_f0;
Mm0=M0*nu_m0;

k_sigma_f0 = 0.1; 
k_sigma_m0 = 0.1;

time = days;       % dt=1/n_dt time unit
n_time = n_dt*time+1;
init_dmg_t = 10; % G&R simulation w/o damage until init_dmg_t
dt =1/n_dt;     % time interval
num_DL= floor(age_max/dt)+1; % number of data per life span

% parameters in nondimensionalization

A_11=1.0/(P_a*r_h);  %See non-dimensinalization in the note
beta = 1/P_a;        %  normalized pressure = beta * Pressure


%iteration parameter
max_itr = 500;
max_more_it=3;        %iteration for rate of mass production
tol =1.0*10^(-7);

z_0 = length;
l_elm = z_0/n_elm; %Length of elements


switch(type)
case 'linear'
        n_node = n_elm+1;
        sh = l_elm;
        npe = 2;
        Data = zeros(n_node, 2);
        for i=1:n_node
            Data(i,1)= sh*(i-1);
            Data(i,2)= 1.0; %nornalized radius   r/r_h
        end
case 'quadratic'
       n_node = 2*n_elm +1;
       sh= l_elm/2.0;
       npe = 3;
       Data = zeros(n_node,2);
       for i=1:n_elm
           ii=2*i-1;
           
           Data(ii,1)=sh*(ii-1);
           Data(ii,2)=1.0;  %normalized radius

           Data(ii+1,1)= sh*ii;
           Data(ii+1,2)=1.0;
       end
       Data(n_node,1)= sh*(n_node-1);
       Data(n_node,2)=1.0;
otherwise
    error('Unknown type of shape function');
end

% Gauss points and weight values for integration
Gauss_Pt=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
   -0.5773502692, 0.5773502692, 0.0, 0.0, 0.0, 0.0;
   -0.774596669241, 0.0, 0.774596669241, 0.0, 0.0, 0.0;
   -0.861136311594, -0.339981043584, 0.339981043584,0.861136311594, 0.0, 0.0;
   -0.906179845938, -0.538469310105, 0.0, 0.538469310105, 0.906179845938, 0.0;  
   -0.9324695142, -0.6612093865, -0.2386191861,0.2386191861, 0.6612093865, ...
       0.9324695142];
Gauss_Wi=[2.000, 0.0, 0.0, 0.0, 0.0, 0.0;
    1.000, 1.000, 0.0, 0.0, 0.0, 0.0;
    0.555555555555, 0.888888888889, 0.555555555555, 0.0, 0.0, 0.0;
    0.347854845137, 0.652145154862, 0.652145154862, 0.347854845137, 0.0, 0.0;
    0.236926885056, 0.478628670499, 0.568888888889,0.478628670499,...
        0.236926885056, 0.0;
    0.1713244924, 0.3607615730, 0.4679139346,0.4679139346, 0.3607615730,...
        0.1713244924];

num_pts=n_elm*n_Gauss_pt ;
TD_L1=zeros(num_pts,n_time+1);    % lambda_1(points,time)
TD_L2=zeros(num_pts,n_time+1);    % lambda_2(points, time)

TD_mc1=zeros(num_pts, n_time+1);   % mass production of collagen fiber in circumferential direction
TD_mc2=zeros(num_pts, n_time+1);   % mass production of collagen fiber in axial direction
TD_mc3=zeros(num_pts, n_time+1);   % mass production of collagen fiber in diagonal direction
TD_mm=zeros(num_pts, n_time+1);   % mass production of SM

%TD_ang=zeros(num_pts, n_time+1);   % mass production of collagen fiber in
%theta
%TD_Zc=zeros(1,n_time+1);    %thickness  at the given point (Z=Zc)

% initial conditions
TD_L1(:,1)=1.0;
TD_L2(:,1)=1.0;

%survival fraction Q(t)
DQ_c = zeros(num_DL,1);
DQ_m = zeros(num_DL,1);
DQ_c(1)=0.0;
DQ_m(1)=0.0;
for i=2:num_DL
    t= dt*(i-1);
    DQ_c(i) = DQ_c(i-1) + 0.5*dt*(q_i(0, t-dt)+q_i(0,t));
    DQ_m(i) = DQ_m(i-1) + 0.5*dt*(q_i(1, t-dt)+q_i(1,t));
end

mean_age_m = DQ_m(num_DL);
mean_age_c = DQ_c(num_DL);

for i=1:num_DL
    DQ_m(i) = (mean_age_m-DQ_m(i))/mean_age_m;
    DQ_c(i) = (mean_age_c-DQ_c(i))/mean_age_c;
end

%initial mass production m^k(t=0+), 
%initial fiber alingment  alpha_0 (t=0+) 

cnt=0;
m_basal_k = Mf0/mean_age_c;
m_basal_m = Mm0/mean_age_m;

for elm=1:n_elm
         for k=1:n_Gauss_pt
             cnt=cnt+1;

             TD_mc1(cnt,1)=m_basal_k(1);
             TD_mc2(cnt,1)=m_basal_k(2);
             TD_mc3(cnt,1)=m_basal_k(3);
             TD_mm(cnt,1)=m_basal_m;
        end
end

n_Zc=floor(num_pts/2.0);                                  %  point in the middle of z-axis

%initial guess
x_n=zeros(n_node*2,1);     % (n)th values for (r,z)
r_nod=zeros(n_node,1);
z_nod=zeros(n_node,1);
for i=1:n_node
    x_n(i*2-1)= Data(i,2);       %guess initial radius
    x_n(i*2)=Data(i,1);          %guess initial z
end

x_np=zeros(n_node*2,1);   % (n-1)th values of r and z
x_pr1 = zeros(n_node*2,1); % previous time value
x_pr2 = zeros(n_node*2,1);

for i=1:n_node*2          %initial guess
     x_pr1(i) = x_n(i);
     x_pr2(i) = x_pr1(i);
end

glk = zeros(n_node*2, n_node*2);
glf = zeros(n_node*2,1);

itnum=0;
p_count=0;

%--------------------------------------------------
% Begining of the loops
%--------------------------------------------------

P=P_a;
p_count=0;

for nt=1:n_time
   
 if nt==1
    TD_mc1(:,nt+1) = TD_mc1(:,nt);   %initial guess for the mass
    TD_mc2(:,nt+1) = TD_mc2(:,nt);
    TD_mc3(:,nt+1) = TD_mc3(:,nt);
    TD_mm(:, nt+1) = TD_mm(:, nt);
 else
    TD_mc1(:,nt+1) = 2.0*TD_mc1(:,nt)-TD_mc1(:, nt-1);   %initial guess for the mass
    TD_mc2(:,nt+1) = 2.0*TD_mc2(:,nt)-TD_mc2(:, nt-1);
    TD_mc3(:,nt+1) = 2.0*TD_mc3(:,nt)-TD_mc3(:, nt-1);
    TD_mm(:,nt+1) = 2.0*TD_mm(:,nt)-TD_mm(:, nt-1);
 end
 for i=1:n_node*2
     x_n(i) = 2.0*x_pr1(i)-x_pr2(i);          %initial guess
 end
 % guess for the mass production rate, then find stretch
 % ,stresses, and total mass, then new guess for the mass production
 % iterate the loop 'max_more_it' times
 
 current_t=dt*nt;
 message = sprintf('%3.1f days, pressure %3.4f mmHg \n', current_t, P/133.322668);
 disp(message);
 
 
 for more_it=1:max_more_it   % interation for mass
     
 if itnum>max_itr
      error('Too many itteration! teminate the program');
      pause
 end
 
%itnum  %print itnum on the console
 
 itnum=0;    %number of iteration
 itr_error = 10.0;  %initial value of iteration error

 while ( itnum <= max_itr & itr_error > tol)
     itnum = itnum+1;
     
     if (current_t >= op_time) && (p_count<20)
         p_count=p_count+1;
         P = (1-p_count/20)*P_a+(p_count/20)*0.4*P_a;
         k_stent = 10*p_count;
     end
                    
     %initializing the matrix
     for i=1:n_node*2
         glf(i)=0.0;
         x_np(i)=x_n(i);       %store phi(n) value to phi(n-1)
         x_n(i)=0.0;              % and initialize phi(n)
         for j=1:n_node*2
              glk(i, j)=0.0;
         end
     end
             
     %integration using Gauss-integral
     count=0;
     for elm=1:n_elm
         %evaluating phi and phi' at a gauss point on every elements
         for k=1:n_Gauss_pt
             count=count+1;
             Z=((elm-1)+(Gauss_Pt(n_Gauss_pt,k)+1)/2.0)*l_elm;
             
             R=0.0;
             dR=0.0;
             
             z= 0.0;
             r=0.0;
             dz =0.0;
             dr = 0.0;
           
             for i=1:npe
                     index=(elm-1)*(npe-1)+i;
                     R = R+Data(index,2)*shp1d(index,0,Z);
                     dR = dR+Data(index,2)*shp1d(index,1,Z);
                     
                     z = z +x_np(index*2)*shp1d(index,0,Z);
                     r = r +x_np(index*2-1)*shp1d(index,0,Z);
                     dz =dz+x_np(index*2)*shp1d(index,1,Z);
                     dr = dr+x_np(index*2-1)*shp1d(index,1,Z);
                     
             end

            sz2r2=sqrt(dz^2+dr^2);
            s1R2=sqrt(1+dR^2);
            
            L1=r/R;                                                %lambda_1 (circumflex)
            L2=sz2r2/s1R2;     %lambda_theta (z-direction)
            TD_L1(count, nt+1) = L1;
            TD_L2(count, nt+1) = L2;
               
            dL1dr=1/R;
            dL2dr2=dr/(s1R2*sz2r2);
            dL2dz2=dz/(s1R2*sz2r2);
            ddL2ddr2=dz^2 / (s1R2*sz2r2^3);
            ddL2ddz2=dr^2 / (s1R2*sz2r2^3);
            ddL2dr2dz2 = - dr*dz / (s1R2*sz2r2^3);
                
            if nt+1 <= num_DL  %If there still remains initial mass 
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
     
                Mf_r=Mf0*DQ_c(nt+1); 
                Mt = Mf_r(1)+Mf_r(2)+Mf_r(3)+Mf_r(4);
                
                dwdL_k1= A_11*Mf_r(1)*G_h*dWcdx(Ln_k1);
                dwdL_k2= A_11*Mf_r(2)*G_h*dWcdx(Ln_k2);
                dwdL_k3= A_11*Mf_r(3)*G_h*dWcdx(Ln_k3);
                dwdL_k4= A_11*Mf_r(4)*G_h*dWcdx(Ln_k4);
                
                
                ddwddL_k1=A_11*Mf_r(1)*G_h^2*ddWcddx(Ln_k1);
                ddwddL_k2=A_11*Mf_r(2)*G_h^2*ddWcddx(Ln_k2);
                ddwddL_k3=A_11*Mf_r(3)*G_h^2*ddWcddx(Ln_k3);
                ddwddL_k4=A_11*Mf_r(4)*G_h^2*ddWcddx(Ln_k4);
                
                dwdL1=dwdL_k1+dwdL_k3*dL_k3dL1+dwdL_k4*dL_k4dL1;
                dwdL2=dwdL_k2+dwdL_k3*dL_k3dL2+dwdL_k4*dL_k4dL2;
                ddwddL1=ddwddL_k1+ddwddL_k3*dL_k3dL1^2+dwdL_k3*ddL_k3ddL1 + ...
                    ddwddL_k4*dL_k4dL1^2+dwdL_k4*ddL_k4ddL1;
                ddwddL2=ddwddL_k2+ddwddL_k3*dL_k3dL2^2+dwdL_k3*ddL_k3ddL2 + ...
                    ddwddL_k4*dL_k4dL2^2+dwdL_k4*ddL_k4ddL2;
                ddwdL1dL2=ddwddL_k3*dL_k3dL1*dL_k3dL2+dwdL_k3*ddL_k3dL1dL2+ ...
                    ddwddL_k4*dL_k4dL1*dL_k4dL2+dwdL_k4*ddL_k4dL1dL2;  
                
                % we assume that the initial SM removed with the elastin
                % degradation
                
                Mm = Mm0*DQ_m(nt+1)*Rm_exp(Z, z_0, current_t, init_dmg_t);
                Mt = Mt+Mm;
                Ln_m = G_m * L1;
                dwdL1 = dwdL1+ A_11*Mm*G_m*dWmdx(Ln_m);
                ddwddL1 = ddwddL1 + A_11*Mm*G_m^2*ddWmddx(Ln_m);
                
            else  % All of initial fiber families and SM are removed 
                dwdL1=0.0;
                dwdL2=0.0;
                ddwddL1=0.0;
                ddwddL2=0.0;
                ddwdL1dL2=0.0;
                Mt = 0.0;
                Mm = 0.0;
            end
            
            % strain energy due to elastin layer
            
            Me = Me0*Rm_exp(Z, z_0, current_t, init_dmg_t);
            Mt = Mt+Me;
            
            Ln_e1 = G_e(1)*L1;
            Ln_e2 = G_e(2)*L2;
                
            dwdL_e1= A_11*Me*G_e(1)*dWedx(1,Ln_e1, Ln_e2);
            dwdL_e2= A_11*Me*G_e(2)*dWedx(2,Ln_e1, Ln_e2);
            ddwddL_e1=A_11*Me*(G_e(1))^2*ddWeddx(1,Ln_e1,Ln_e2);
            ddwddL_e2=A_11*Me*(G_e(2))^2*ddWeddx(2,Ln_e1,Ln_e2);
            ddwddL_e1dL_e2 = A_11*Me*G_e(1)*G_e(2)*ddWeddx(3,Ln_e1,Ln_e2);
  
            dwdL1 = dwdL1 +dwdL_e1;
            dwdL2 = dwdL2 +dwdL_e2;
            ddwddL1 = ddwddL1 +ddwddL_e1;
            ddwddL2 = ddwddL2 +ddwddL_e2;
            ddwdL1dL2 = ddwdL1dL2 +ddwddL_e1dL_e2;
            
            %-----------------------------------------------------------
            %       numerical integration
            % int(0 to t) [ ff(tau)]dtau = ff(0)*0.5*dtau+
            % ff(dtau)*dtau+ff(2*dtau)*dtau+ . .+ff(t)*0.5*dtau
            %-----------------------------------------------------------
            if nt+1 <= num_DL
                n_tau0 = 1;
            else
                n_tau0= nt-num_DL+2;
            end
            for n_tau=n_tau0: (nt+1) 
                   tau=(n_tau-1)*dt;
                   L1_tau=TD_L1(count,n_tau);
                   L2_tau=TD_L2(count,n_tau);
                   
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
                           TD_mc3(count, n_tau), TD_mc3(count, n_tau)]*q_i(0, current_t-tau);
                   mm_tau = TD_mm(count, n_tau) * q_i(1, current_t-tau)...
                       *Rm_exp(Z, z_0, current_t, init_dmg_t);  
                   Mm = Mm+mm_tau*w_dt;
                   Mt = Mt+ (mc_tau(1)+mc_tau(2)+mc_tau(3)+mc_tau(4))*w_dt...
                       +mm_tau*w_dt;
                   
                   Ln_k1 = G_h*L1 / L1_tau;
                   Ln_k2 = G_h*L2 / L2_tau;
                   Ln_k3 = G_h*L_k3 / L_k3_tau;
                   Ln_k4 = G_h*L_k4 / L_k4_tau;
                   Ln_m =  G_m*L1/L1_tau;
            
                   dwdL_k1 = A_11* mc_tau(1) * (G_h/L1_tau)*dWcdx(Ln_k1);
                   dwdL_k2 = A_11* mc_tau(2) * (G_h/L2_tau)*dWcdx(Ln_k2);
                   dwdL_k3 = A_11* mc_tau(3) * (G_h/L_k3_tau)*dWcdx(Ln_k3);
                   dwdL_k4 = A_11* mc_tau(4) * (G_h/L_k4_tau)*dWcdx(Ln_k4);                 
                   
                   ddwddL_k1 = A_11*mc_tau(1) * (G_h/L1_tau)^2*ddWcddx(Ln_k1);
                   ddwddL_k2 = A_11*mc_tau(2) * (G_h/L2_tau)^2*ddWcddx(Ln_k2);
                   ddwddL_k3 = A_11*mc_tau(3) * (G_h/L_k3_tau)^2*ddWcddx(Ln_k3);
                   ddwddL_k4 = A_11*mc_tau(4) * (G_h/L_k4_tau)^2*ddWcddx(Ln_k4);
                   
                   dwdL1=dwdL1 + (dwdL_k1+dwdL_k3*dL_k3dL1+dwdL_k4*dL_k4dL1)*w_dt; 
                   dwdL2=dwdL2 + (dwdL_k2+dwdL_k3*dL_k3dL2+dwdL_k4*dL_k4dL2)*w_dt;
                   ddwddL1=ddwddL1+ (ddwddL_k1+ddwddL_k3*dL_k3dL1^2+dwdL_k3*ddL_k3ddL1 + ...
                          ddwddL_k4*dL_k4dL1^2+dwdL_k4*ddL_k4ddL1)*w_dt ;
                   ddwddL2=ddwddL2+ (ddwddL_k2+ddwddL_k3*dL_k3dL2^2+dwdL_k3*ddL_k3ddL2 + ...
                          ddwddL_k4*dL_k4dL2^2+dwdL_k4*ddL_k4ddL2)*w_dt ;
                   ddwdL1dL2=ddwdL1dL2+(ddwddL_k3*dL_k3dL1*dL_k3dL2+dwdL_k3*ddL_k3dL1dL2+ ...
                          ddwddL_k4*dL_k4dL1*dL_k4dL2+dwdL_k4*ddL_k4dL1dL2) * w_dt; 
                      
                   %Contribution of SM
                   dwdL1 = dwdL1+A_11* mm_tau * (G_m/L1_tau)*dWmdx(Ln_m)*w_dt;
                   ddwddL1 = ddwddL1+A_11*mm_tau * (G_m/L1_tau)^2*ddWmddx(Ln_m)*w_dt;
            end
  
            f_act = 1-((La_M-1.0)/(La_M-La_0))^2;   
            dwdL1_act = A_11*Mm*(1.0/L1)*(Sa/rho)*f_act;
            %ddwddL1_act = A_11*Mm*Sa/rho*(1/L1_s)^2*2.0*(La_M-L1/L1_s)/(La_M-La_0)^2;
            ddwddL1_act = -A_11*Mm*(1.0/L1^2)*(Sa/rho)*f_act;
            
            dwdL1 = dwdL1+dwdL1_act;
            ddwddL1 = ddwddL1+ddwddL1_act;
            
            dwdr =dwdL1*dL1dr;
            dwdr2=dwdL2*dL2dr2;
            dwdz2=dwdL2*dL2dz2;
            ddwddr = ddwddL1*dL1dr^2;
            ddwddr2= ddwddL2*dL2dr2^2+dwdL2* ddL2ddr2;
            ddwddz2= ddwddL2*dL2dz2^2+dwdL2*ddL2ddz2;
            ddwdrdr2=ddwdL1dL2*dL1dr*dL2dr2;
            ddwdrdz2=ddwdL1dL2*dL1dr*dL2dz2;
            ddwdr2dz2 = ddwddL2*dL2dr2*dL2dz2+dwdL2*ddL2dr2dz2;

            %Global K , K_ij
            WI=Gauss_Wi(n_Gauss_pt, k)*l_elm/2;
            for i=1:npe
                ii = (elm-1)*(npe-1)+i;
                ii_1=ii*2-1;
                ii_2=ii*2;
                Si = shp1d(ii,0,Z);
                dSi=shp1d(ii,1,Z);
                
                
                glf(ii_1)=glf(ii_1)+2*pi*((dwdr*Si + dwdr2*dSi)*R*s1R2- ...
                    beta*P*dz*r*Si)*WI;
                
                if current_t>=op_time
                    if r_stent > r
                        F_penalty = -1.0*sz2r2*k_stent*(r_stent-r)*Si;
                        glf(ii_1)=glf(ii_1)+2*pi*F_penalty*WI;
                    end
                end
                glf(ii_2)=glf(ii_2)+2*pi*(dwdz2*dSi*R*s1R2 + ...
                    beta*P*dr*r*Si)*WI;
                
                for j=1:npe
                    jj= (elm-1)*(npe-1)+j;
                    jj_1=jj*2-1;
                    jj_2=jj*2;
                    SiSj=shp1d(ii,0,Z)*shp1d(jj,0,Z);
                    SidSj=shp1d(ii,0,Z)*shp1d(jj,1,Z);
                    dSiSj=shp1d(ii,1,Z)*shp1d(jj,0,Z);
                    dSidSj=shp1d(ii,1,Z)*shp1d(jj,1,Z);
                    
                    
                    glk(ii_1,jj_1)=glk(ii_1,jj_1)+2*pi*((ddwddr*SiSj+ddwdrdr2*(SidSj+dSiSj)...
                        +ddwddr2*dSidSj)*R*s1R2-beta*P*dz*SiSj)*WI;
                    glk(ii_1,jj_2)=glk(ii_1,jj_2)+2*pi*((ddwdrdz2*SidSj+ ddwdr2dz2*dSidSj)*...
                        R*s1R2-beta*P*r*SidSj)*WI;
                    
                    if current_t>=op_time
                        if r_stent > r
                            K11_penalty = (sz2r2*k_stent*SiSj-(dr/sz2r2)*k_stent*(r_stent-r)*SidSj);
                            K12_penalty =-1*(dz/sz2r2)*k_stent*(r_stent-r)*SidSj;
                            glk(ii_1,jj_1)=glk(ii_1,jj_1)+2*pi*K11_penalty*WI;
                            glk(ii_1,jj_2)=glk(ii_1,jj_2)+2*pi*K12_penalty*WI;
                        end
                    end
                    glk(ii_2,jj_2)=glk(ii_2,jj_2)+2*pi*(ddwddz2*dSidSj*R*s1R2+beta*P*r*SidSj)*WI;
                    
                    glk(ii_2,jj_1)=glk(ii_2,jj_1)+2*pi*((ddwdrdz2*dSiSj+ddwdr2dz2*dSidSj)*...
                        R*s1R2+beta*P*(dr*SiSj+r*SidSj))*WI;
                    
                end
            end       
        end
    end
       
    %Imposing boundary conditions
    for i=1:n_node*2
        glk(2,i)=0.0;
        glk(n_node*2,i)=0.0;
    end
    glf(2) = 0.0;
    glk(2,2)=1.0;
    glf(n_node*2)=0;
    glk(n_node*2, n_node*2)=1.0;
  
    if bc==1
        for i=1:n_node*2
            glk(1,i)=0.0;
            glk(n_node*2-1,i)=0.0;
        end
        glf(1) = 0.0;
        glk(1,1)=1.0;
        glf(n_node*2-1)=0;
        glk(n_node*2-1, n_node*2-1)=1.0;
    end
    %glf
    %glk
    %solving matrix
    
    dx=glk\glf;
    
    x_n=x_np-dx;
    % error
    itr_error = sqrt(sum(dx.^2)/sum(x_n.^2))
 end
 
 %the end of FEM calculation for each pressure P
 % Save the results
 
 % using the FEM solution, 
 % calculate stress, strain, etc. at the Gauss point
 
 %print data into a file at every 'print_step'
 
 if (more_it==max_more_it && (mod(nt,print_step)==1 || nt==init_dmg_t*n_dt-1))
     p_count=p_count + 1;
     string1=name;
     print_day = floor((nt-1)/n_dt);
     add_str = sprintf('_s%i.dat',print_day);       
     outfile_name =strcat(string1,add_str);
     fid_out_1=fopen(outfile_name,'w');
     out_string = sprintf('Time=%f\nWriting ''%s''',print_day, outfile_name);
     disp(out_string)
 end
 
 count=0;
 

 
 for elm=1:n_elm
     for k=1:n_Gauss_pt
         count=count+1;
         Z=((elm-1)+(Gauss_Pt(n_Gauss_pt,k)+1)/2.0)*l_elm;
         
         R=0.0;
         dR=0.0;
         
         z= 0.0;
         r=0.0;
         dz =0.0;
         dr = 0.0;
         
         for i=1:npe
                 index=(elm-1)*(npe-1)+i;
                 R = R+Data(index,2)*shp1d(index,0,Z);
                 dR = dR+Data(index,2)*shp1d(index,1,Z);
                 
                 z = z +x_np(index*2)*shp1d(index,0,Z);
                 r = r +x_np(index*2-1)*shp1d(index,0,Z);
                 dz =dz+x_np(index*2)*shp1d(index,1,Z);
                 dr = dr+x_np(index*2-1)*shp1d(index,1,Z);
                
         end

        sz2r2=sqrt(dz^2+dr^2);
        s1R2=sqrt(1+dR^2);
            
        L1=r/R;                   %lambda_1 (circumflex)
        L2=sz2r2/s1R2;            %lambda_theta (z-direction)
        
        TD_L1(count, nt+1) = L1;
        TD_L2(count, nt+1) = L2;
               
        dL1dr=1/R;
        dL2dr2=dr/(s1R2*sz2r2);
        dL2dz2=dz/(s1R2*sz2r2);
        ddL2ddr2=dz^2 / (s1R2*sz2r2^3);
        ddL2ddz2=dr^2 / (s1R2*sz2r2^3);
        ddL2dr2dz2 = - dr*dz / (s1R2*sz2r2^3);
           
         if nt+1 <= num_DL  %If there is still initial mass 
                L_k3=sqrt( (L1*sin(phi0))^2+(L2*cos(phi0))^2);
                L_k4=sqrt( (L1*sin(pi-phi0))^2+(L2*cos(pi-phi0))^2);
            
                dL_k3dL1=L1*sin(phi0)^2 / L_k3;  %dL_k4dL1=dL_k3dL1
                dL_k3dL2=L2*cos(phi0)^2 / L_k3;
                dL_k4dL1=L1*sin(pi-phi0)^2 / L_k4;
                dL_k4dL2=L2*cos(pi-phi0)^2 / L_k4;
                
                Ln_k1 = G_h*L1;
                Ln_k2 = G_h*L2;
                Ln_k3 = G_h*L_k3;
                Ln_k4 = G_h*L_k4;
     
                Mf_r=Mf0*DQ_c(nt+1); 
                Mf = Mf_r(1)+Mf_r(2)+Mf_r(3)+Mf_r(4);
                
                dwdL_k1= Mf_r(1)*G_h*dWcdx(Ln_k1);
                dwdL_k2= Mf_r(2)*G_h*dWcdx(Ln_k2);
                dwdL_k3= Mf_r(3)*G_h*dWcdx(Ln_k3);
                dwdL_k4= Mf_r(4)*G_h*dWcdx(Ln_k4);
                
                dwdL1_f=dwdL_k1+dwdL_k3*dL_k3dL1+dwdL_k4*dL_k4dL1;
                dwdL2_f=dwdL_k2+dwdL_k3*dL_k3dL2+dwdL_k4*dL_k4dL2;
                
                Mm = Mm0*DQ_m(nt+1)*Rm_exp(Z, z_0, current_t, init_dmg_t);
                
                Ln_m = G_m * L1;
                dwdL1_m =  Mm*G_m*dWmdx(Ln_m);
        else
                dwdL1_f = 0.0;
                dwdL2_f = 0.0;
                dwdL1_m = 0.0;
                Mf = 0.0;
                Mm = 0.0;
        end
        
        % stress due to elastin layer
            
         Me = Me0*Rm_exp(Z, z_0, current_t, init_dmg_t);
         
         Ln_e1 = G_e(1)*L1;
         Ln_e2 = G_e(2)*L2;
                
         dwdL_e1= Me*G_e(1)*dWedx(1,Ln_e1, Ln_e2);
         dwdL_e2= Me*G_e(2)*dWedx(2,Ln_e1, Ln_e2);
         
         
     
       %numerical integration
        if nt+1 <= num_DL
            n_tau0 = 1;
        else
            n_tau0= nt-num_DL+2;
        end
        for n_tau=n_tau0:nt+1  
               tau=(n_tau-1)*dt;
               L1_tau=TD_L1(count,n_tau);
               L2_tau=TD_L2(count,n_tau);
                   
               phi = phi0;  % phi(\tau) is not changing through time
                   
               L_k3=sqrt( (L1*sin(phi))^2+(L2*cos(phi))^2);
               L_k4=sqrt( (L1*sin(pi-phi))^2+(L2*cos(pi-phi))^2);
            
               dL_k3dL1=L1*sin(phi)^2 / L_k3;  %dL_k4dL1=dL_k3dL1
               dL_k3dL2=L2*cos(phi)^2 / L_k3;
               dL_k4dL1=L1*sin(pi-phi)^2 / L_k4;
               dL_k4dL2=L2*cos(pi-phi)^2 / L_k4;
                   
               L_k3_tau=sqrt( (L1_tau*sin(phi))^2+(L2_tau*cos(phi))^2);
               L_k4_tau = L_k3_tau;
                   
               if ((n_tau==n_tau0) || (n_tau==nt+1))
                   w_dt=dt*0.5; 
               else
                   w_dt=dt;
               end
                   
               mc_tau=[TD_mc1(count,n_tau), TD_mc2(count, n_tau), ...
                       TD_mc3(count, n_tau), TD_mc3(count, n_tau)]*q_i(0, current_t-tau);
               mm_tau = TD_mm(count, n_tau) * q_i(1, current_t-tau)...
                   *Rm_exp(Z, z_0, current_t, init_dmg_t);;  
               Mm = Mm+mm_tau*w_dt;
               Mf = Mf+ (mc_tau(1)+mc_tau(2)+mc_tau(3)+mc_tau(4))*w_dt;
                   
               Ln_k1 = G_h*L1 / L1_tau;
               Ln_k2 = G_h*L2 / L2_tau;
               Ln_k3 = G_h*L_k3 / L_k3_tau;
               Ln_k4 = G_h*L_k4 / L_k4_tau;
            
               dwdL_k1 = mc_tau(1) * (G_h/L1_tau)*dWcdx(Ln_k1);
               dwdL_k2 = mc_tau(2) * (G_h/L2_tau)*dWcdx(Ln_k2);
               dwdL_k3 = mc_tau(3) * (G_h/L_k3_tau)*dWcdx(Ln_k3);
               dwdL_k4 = mc_tau(4) * (G_h/L_k4_tau)*dWcdx(Ln_k4);
                   
               dwdL1_f = dwdL1_f + (dwdL_k1+dwdL_k3*dL_k3dL1+dwdL_k4*dL_k4dL1)*w_dt; 
               dwdL2_f = dwdL2_f + (dwdL_k2+dwdL_k3*dL_k3dL2+dwdL_k4*dL_k4dL2)*w_dt; 
               
               %Contribution of SM
               Ln_m = G_m*L1/L1_tau;
               dwdL1_m = dwdL1_m + mm_tau * (G_m/L1_tau)*dWmdx(Ln_m)*w_dt;
        end 
             
        dwdL1 = dwdL1_f +dwdL_e1+dwdL1_m;
        dwdL2 = dwdL2_f +dwdL_e2;
          
        Mt = Mf+Me+Mm;
        
        T_f_p1 = (1.0/L2)*dwdL1_f;
        T_f_p2 = (1.0/L1)*dwdL2_f;
        
        T_p1 = (1.0/L2)*dwdL1;
        T_p2 = (1.0/L1)*dwdL2;
        
        %f_act = 1-((La_M-L1)/(La_M-La_0))^2;
        f_act = 1-((La_M-1.0)/(La_M-La_0))^2;
        dwdL1_act = (1.0/L1)*Mm*(Sa/rho)*f_act; 
        
        T_m = (1.0/L2)*(dwdL1_m + dwdL1_act);
            
        T_1 =  T_p1+(1.0/L2)*dwdL1_act;
        T_2 =  T_p2;
 
        TD_L1(count, nt+1) = L1;
        TD_L2(count, nt+1) = L2;
        
        thick_pt = Mt /(rho*L1*L2);
        thick_f_pt = Mf/(rho*L1*L2);
        thick_m_pt = Mm/(rho*L1*L2);
        
        pS1= T_f_p1 / thick_f_pt;
        pS2= T_f_p2 / thick_f_pt;
        if thick_m_pt > 0
            pSm= T_m / thick_m_pt;
        else
            pSm =0.0;
        end
        
        S1 =  T_1/thick_pt;
        S2 =  T_2/thick_pt;
        
        phi = atan( L1*sin(phi0)/(L2*cos(phi0)));
        
        sigma_1=pS1;
        sigma_2=pS2;
        sigma_3=sqrt((pS1*sin(phi))^2+(pS2*cos(phi))^2);
        sigma_m = pSm;
        
        if nt<init_dmg_t
            k_f=k_sigma_f0;
            k_m=k_sigma_m0;
        else
            k_f=k_sigma_f;
            k_m=k_sigma_m;
        end
        
        % m_i(n_c, m_basal,dn_stress, Kc)
        Mf0_total=sum(Mf0);
        
        TD_mc1(count,nt+1)=m_i(Mf/Mf0_total,m_basal_k(1), ...
            (sigma_1-sigma_f0)/sigma_f0, k_f);    %sigma/sigma_0 = normalized value
        TD_mc2(count,nt+1)=m_i(Mf/Mf0_total,m_basal_k(2), ...
            (sigma_2-sigma_f0)/sigma_f0, k_f);
        TD_mc3(count,nt+1)=m_i(Mf/Mf0_total,m_basal_k(3), ...
            (sigma_3-sigma_f0)/sigma_f0, k_f);
        TD_mm(count,nt+1)=m_i(Mm/Mm0,m_basal_m,  ...
            (sigma_m-sigma_m0)/sigma_m0,k_m);
        
        %if count== floor(n_elm/10)*n_Gauss_pt 
        %    [L1 L2 L_k3]
        %    [Mf/Mf0_total  Mm/Mm0]
        %   [(sigma_1-sigma_f0)/sigma_f0, (sigma_3-sigma_f0)/sigma_f0, (sigma_m-sigma_m0)/sigma_m0]
        %   [TD_mc1(count,nt+1)/m_basal_k(1), TD_mc3(count,nt+1)/m_basal_k(3), TD_mm(count,nt+1)/m_basal_m]
        %end 
     
        %record values if time is integer 
        if ( more_it==max_more_it && (mod(nt,print_step)== 1 || nt==init_dmg_t*n_dt-1))
            if flag_LinearC > 0
                fprintf(fid_out_1,'%10.9f\t%10.9f\t%12.10f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\n',...
                    z*r_h, r*r_h, thick_pt, S1, S2, Mt, Mm, Me);

            else
                % calculating linearized material parameters
                if nt+1 <= num_DL  %If there is still initial mass
                    Ln_k1 = G_h *L1;
                    Ln_k2 = G_h *L2;
                    Ln_k3 = G_h *L_k3;
                    Ln_k4 = G_h *L_k4;
                    Ln_m = G_m * L1;

                    dW_k1dC_11 = (0.5/Ln_k1)*dWcdx(Ln_k1) * G_h^2;
                    dW_k3dC_11 = (0.5/Ln_k3)*dWcdx(Ln_k3) * G_h^2* sin(phi0)^2;
                    dW_k4dC_11 = (0.5/Ln_k4)*dWcdx(Ln_k4) * G_h^2* sin(pi-phi0)^2;

                    dW_mdC_11 = (0.5/Ln_m)*dWmdx(Ln_m) * G_m^2;

                    dW_k2dC_22 = (0.5/Ln_k2)*dWcdx(Ln_k2) * G_h^2;
                    dW_k3dC_22 = (0.5/Ln_k3)*dWcdx(Ln_k3) * G_h^2* cos(phi0)^2;
                    dW_k4dC_22 = (0.5/Ln_k4)*dWcdx(Ln_k4) * G_h^2* cos(pi-phi0)^2;

                    Mf_r=Mf0*DQ_c(nt+1);
                    Mf = Mf_r(1)+Mf_r(2)+Mf_r(3)*2.0;
                    Mm = Mm0*DQ_m(nt+1)*Rm_exp(Z, z_0,current_t, init_dmg_t);


                    dWdC_11 = (rho/Mt)*(Mf_r(1)*dW_k1dC_11+ ...
                        Mf_r(3)*dW_k3dC_11+Mf_r(4)*dW_k4dC_11 +Mm*dW_mdC_11);
                    dWdC_22 = (rho/Mt)*(Mf_r(2)*dW_k2dC_22+ ...
                        Mf_r(3)*dW_k3dC_22+Mf_r(4)*dW_k4dC_22);
                    dWdC_33 = 0.0;

                    ddW_k1ddLn_k1_2 = (kc(2)/2.0)*(1.0+2.0*kc(3)*(Ln_k1^2-1)^2)*exp(kc(3)*(Ln_k1^2-1)^2);
                    ddW_k2ddLn_k2_2 = (kc(2)/2.0)*(1.0+2.0*kc(3)*(Ln_k2^2-1)^2)*exp(kc(3)*(Ln_k2^2-1)^2);
                    ddW_k3ddLn_k3_2 = (kc(2)/2.0)*(1.0+2.0*kc(3)*(Ln_k3^2-1)^2)*exp(kc(3)*(Ln_k3^2-1)^2);
                    ddW_k4ddLn_k4_2 = (kc(2)/2.0)*(1.0+2.0*kc(3)*(Ln_k4^2-1)^2)*exp(kc(3)*(Ln_k4^2-1)^2);
                    ddW_mddLn_m_2 = (kc(4)/2.0)*(1.0+2.0*kc(5)*(Ln_m^2-1)^2)*exp(kc(5)*(Ln_m^2-1)^2);

                    ddWddC_11 = (rho/Mt)*(Mf_r(1)*ddW_k1ddLn_k1_2*G_h^4 + Mf_r(3)*ddW_k3ddLn_k3_2*(G_h*sin(phi0))^4 ...
                        + Mf_r(4)*ddW_k4ddLn_k4_2*(G_h*sin(pi-phi0))^4+Mm*ddW_mddLn_m_2*G_m^4);
                    ddWddC_22 = (rho/Mt)*(Mf_r(2)*ddW_k2ddLn_k2_2*G_h^4 + Mf_r(3)*ddW_k3ddLn_k3_2*(G_h*cos(phi0))^4 ...
                        + Mf_r(4)*ddW_k4ddLn_k4_2*(G_h*cos(pi-phi0))^4);
                    ddWdC_11dC_22 = (rho/Mt)*(Mf_r(3)*ddW_k3ddLn_k3_2*(G_h*sin(phi0))^2*(G_h*cos(phi0))^2+...
                        Mf_r(4)*ddW_k4ddLn_k4_2*(G_h*sin(pi-phi0))^2*(G_h*cos(pi-phi0))^2);
                    ddWddC_12 = ddWdC_11dC_22;
                    ddWdC_12dC_11 = (rho/Mt)*(Mf_r(3)*ddW_k3ddLn_k3_2*(G_h^2*sin(phi0)*cos(phi0))*(G_h*sin(phi0))^2+...
                        Mf_r(4)*ddW_k4ddLn_k4_2*(G_h^2*sin(pi-phi0)*cos(pi-phi0))*(G_h*sin(pi-phi0))^2);
                    ddWdC_12dC_22 = (rho/Mt)*(Mf_r(3)*ddW_k3ddLn_k3_2*(G_h^2*sin(phi0)*cos(phi0))*(G_h*cos(phi0))^2+...
                        Mf_r(4)*ddW_k4ddLn_k4_2*(G_h^2*sin(pi-phi0)*cos(pi-phi0))*(G_h*cos(pi-phi0))^2);
                else
                    dWdC_11 = 0.0;
                    dWdC_22 = 0.0;
                    dWdC_33 = 0.0;
                    ddWddC_11 = 0.0;
                    ddWddC_22 = 0.0;
                    ddWdC_11dC_22 = 0.0;
                    ddWddC_12 = 0.0;
                    ddWdC_12dC_11 = 0.0;
                    ddWdC_12dC_22 = 0.0;
                    Mf =0.0;
                    Mm =0.0;
                end

                Ln_e1 = G_e(1)*L1;
                Ln_e2 = G_e(2)*L2;

                dWdC_11 = dWdC_11+(rho/Mt)* Me*G_e(1)/(2.0*L1)*dWedx(1, Ln_e1, Ln_e2);
                dWdC_22 = dWdC_22+(rho/Mt)* Me*G_e(2)/(2.0*L2)*dWedx(2, Ln_e1, Ln_e2);


                if nt+1 <= num_DL
                    n_tau0 = 1;
                else
                    n_tau0= nt-num_DL+2;
                end
                for n_tau=n_tau0:nt+1
                    tau=(n_tau-1)*dt;
                    L1_tau=TD_L1(count,n_tau);
                    L2_tau=TD_L2(count,n_tau);

                    phi = phi0;  % phi(\tau) is not changing through time

                    L_k3=sqrt( (L1*sin(phi))^2+(L2*cos(phi))^2);
                    L_k4=sqrt( (L1*sin(pi-phi))^2+(L2*cos(pi-phi))^2);

                    L_k3_tau=sqrt( (L1_tau*sin(phi))^2+(L2_tau*cos(phi))^2);
                    L_k4_tau = L_k3_tau;

                    if ((n_tau==n_tau0) || (n_tau==nt+1))
                        w_dt=dt*0.5;
                    else
                        w_dt=dt;
                    end

                    mc_tau=[TD_mc1(count,n_tau), TD_mc2(count, n_tau), ...
                        TD_mc3(count, n_tau), TD_mc3(count, n_tau)]*q_i(0, current_t-tau);
                    mm_tau = TD_mm(count, n_tau) * q_i(1, current_t-tau)...
                        *Rm_exp(Z, z_0, current_t, init_dmg_t);;
                    Mm = Mm+mm_tau*w_dt;
                    Mf = Mf+ (mc_tau(1)+mc_tau(2)+mc_tau(3)+mc_tau(4))*w_dt;

                    Ln_k1 = G_h*L1 / L1_tau;
                    Ln_k2 = G_h*L2 / L2_tau;
                    Ln_k3 = G_h*L_k3 / L_k3_tau;
                    Ln_k4 = G_h*L_k4 / L_k4_tau;
                    Ln_m = G_m*L1/L1_tau;

                    dW_k1dC_11 = (0.5/Ln_k1)*dWcdx(Ln_k1) * (G_h/L1_tau)^2;
                    dW_k3dC_11 = (0.5/Ln_k3)*dWcdx(Ln_k3) * (G_h/L_k3_tau)^2* sin(phi0)^2;
                    dW_k4dC_11 = (0.5/Ln_k4)*dWcdx(Ln_k4) * (G_h/L_k4_tau)^2* sin(pi-phi0)^2;

                    dW_mdC_11  = (0.5/Ln_m)*dWmdx(Ln_m) * (G_m/L1_tau)^2;

                    dW_k2dC_22 = (0.5/Ln_k2)*dWcdx(Ln_k2) * (G_h/L2_tau)^2;
                    dW_k3dC_22 = (0.5/Ln_k3)*dWcdx(Ln_k3) * (G_h/L_k3_tau)^2* cos(phi)^2;
                    dW_k4dC_22 = (0.5/Ln_k4)*dWcdx(Ln_k4) * (G_h/L_k4_tau)^2* cos(pi-phi)^2;

                    ddW_k1ddLn_k1_2 = (kc(2)/2.0)*(1.0+2.0*kc(3)*(Ln_k1^2-1)^2)*exp(kc(3)*(Ln_k1^2-1)^2);
                    ddW_k2ddLn_k2_2 = (kc(2)/2.0)*(1.0+2.0*kc(3)*(Ln_k2^2-1)^2)*exp(kc(3)*(Ln_k2^2-1)^2);
                    ddW_k3ddLn_k3_2 = (kc(2)/2.0)*(1.0+2.0*kc(3)*(Ln_k3^2-1)^2)*exp(kc(3)*(Ln_k3^2-1)^2);
                    ddW_k4ddLn_k4_2 = (kc(2)/2.0)*(1.0+2.0*kc(3)*(Ln_k4^2-1)^2)*exp(kc(3)*(Ln_k4^2-1)^2);
                    ddW_mddLn_m_2 = (kc(4)/2.0)*(1.0+2.0*kc(5)*(Ln_m^2-1)^2)*exp(kc(5)*(Ln_m^2-1)^2);

                    dWdC_11 = dWdC_11+(rho/Mt)*(mc_tau(1)*dW_k1dC_11+ ...
                        mc_tau(3)*dW_k3dC_11+mc_tau(4)*dW_k4dC_11+mm_tau*dW_mdC_11)*w_dt;
                    dWdC_22 = dWdC_22+(rho/Mt)*(mc_tau(2)*dW_k2dC_22+ ...
                        mc_tau(3)*dW_k3dC_22+mc_tau(4)*dW_k4dC_22)*w_dt;

                    ddWddC_11 = ddWddC_11+(rho/Mt)*(mc_tau(1)*ddW_k1ddLn_k1_2*(G_h/L1_tau)^4 +...
                        mc_tau(3)*ddW_k3ddLn_k3_2*(G_h/L_k3_tau)^4*sin(phi)^4 ...
                        + mc_tau(4)*ddW_k4ddLn_k4_2*(G_h/L_k4_tau)^4*sin(pi-phi)^4 ...
                        + mm_tau*ddW_mddLn_m_2 )*w_dt;
                    ddWddC_22 = ddWddC_22+(rho/Mt)*(mc_tau(2)*ddW_k2ddLn_k2_2*(G_h/L2_tau)^4 + ...
                        mc_tau(3)*ddW_k3ddLn_k3_2*(G_h/L_k3_tau)^4*cos(phi)^4 ...
                        + mc_tau(4)*ddW_k4ddLn_k4_2*(G_h/L_k4_tau)^4*cos(pi-phi)^4)*w_dt;
                    ddWdC_11dC_22 = ddWdC_11dC_22 + ...
                        (rho/Mt)*(mc_tau(3)*ddW_k3ddLn_k3_2*(G_h/L_k3_tau)^4*sin(phi)^2*cos(phi)^2 ...
                        + mc_tau(4)*ddW_k4ddLn_k4_2*(G_h/L_k4_tau)^4*sin(pi-phi)^2*cos(pi-phi)^2)*w_dt;
                    ddWddC_12 = ddWddC_12+ ...
                        (rho/Mt)*(mc_tau(3)*ddW_k3ddLn_k3_2*(G_h/L_k3_tau)^4*sin(phi)^2*cos(phi)^2 ...
                        + mc_tau(4)*ddW_k4ddLn_k4_2*(G_h/L_k4_tau)^4*sin(pi-phi)^2*cos(pi-phi)^2)*w_dt;;
                    ddWdC_12dC_11 = ddWdC_12dC_11+...
                        (rho/Mt)*(mc_tau(3)*ddW_k3ddLn_k3_2*(G_h/L_k3_tau)^4*sin(phi)^3*cos(phi)+...
                        mc_tau(4)*ddW_k4ddLn_k4_2*(G_h/L_k4_tau)^4*sin(pi-phi)^3*cos(pi-phi))*w_dt;
                    ddWdC_12dC_22 = ddWdC_12dC_22+ ...
                        (rho/Mt)*(mc_tau(3)*ddW_k3ddLn_k3_2*(G_h/L_k3_tau)^4*sin(phi)*cos(phi)^3+...
                        mc_tau(4)*ddW_k4ddLn_k4_2*(G_h/L_k4_tau)^4*sin(pi-phi)*cos(pi-phi)^3)*w_dt;
                end

                Mt=Mf+Me+Mm;

                f_act = 1-((La_M-1.0)/(La_M-La_0))^2;
                df_actdL1 = 2*(La_M-1.0)/(La_M-La_0)^2*(R/r);
                dWdC_11_act =(Mm/Mt)*Sa/(2.0*L1)*f_act*(R/r);
                ddWddC_11_act = (Mm/Mt)*Sa/(4.0*L1^2)*(df_actdL1-f_act/L1)*(R/r);

                dWdC_11 = dWdC_11+dWdC_11_act;
                ddWddC_11 = ddWddC_11+ddWddC_11_act;

                T0_11= 2.0*L1^2*dWdC_11;
                T0_22= 2.0*L2^2*dWdC_22;
                T0_33= 2.0/(L1*L2)^2*dWdC_33;

                T0_12=0.0;

                C_1111 = T0_11+4.0*L1^4*ddWddC_11;
                C_2222 = T0_22+4.0*L2^4*ddWddC_22;
                C_1122 = 4*L1^2*L2^2*ddWdC_11dC_22;
                C_2211 = C_1122;
                C_1212 = T0_22+4.0*L1^2*L2^2*ddWddC_12;
                C_2121 = T0_11+4.0*L1^2*L2^2*ddWddC_12;
                C_1221 = 4.0*L1^2*L2^2*ddWddC_12;
                C_2112 = 4.0*L1^2*L2^2*ddWddC_12;
                C_2323 = T0_33;
                C_3232 = T0_22;
                C_3131 = T0_11;
                C_1313 = T0_33;
                % when the fiber oriented symmetrically, C_1211
                % C_1112, C_1222, and C_2212 are zero
                C_1211 = T0_12+4*L1^3*L2*ddWdC_12dC_11;
                C_1112 = T0_12+4*L1^3*L2*ddWdC_12dC_11;
                C_2111 = 4*L1^3*L2*ddWdC_12dC_11;
                C_1121 = 4*L1^3*L2*ddWdC_12dC_11;
                C_2122 = T0_12 + 4*L1*L2^3*ddWdC_12dC_22;
                C_2221 = T0_12 + 4*L1*L2^3*ddWdC_12dC_22;
                C_1222 = 4*L1*L2^3*ddWdC_12dC_22;
                C_2212 = 4*L1*L2^3*ddWdC_12dC_22;


                fprintf(fid_out_1,'%10.9f\t%10.9f\t%12.10f\t%10.5f\t%10.5f\t%10.5f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',...
                    z*r_h, r*r_h, thick_pt, T0_11-T0_33, T0_22-T0_33, T0_12,C_1111, C_2222, C_1122, ...
                    C_2211, C_1212, C_2121, C_1221, C_2112, C_1211, C_1112, ...
                    C_2122, C_2221, C_2323, C_3232, C_3131, C_1313);
            end
        end
     end 
 end
 if ( more_it==max_more_it && (mod(nt,print_step)== 1 || nt==init_dmg_t*n_dt-1))
     fclose(fid_out_1);
 end
end


 for i=1:n_node*2
     x_pr2(i) = x_pr1(i);
     x_pr1(i) = x_n(i);                      %initial guess
 end
 
 if mod(floor(current_t-dt),100)==0
     for i=1:n_node
         r_nod(i) = x_n(i*2-1);
         z_nod(i) = x_n(i*2);
     end
     plot(z_nod, r_nod); hold on;
 end
 
end

save(name)


clear *
