% AAA_main.m : the main file for G&R of an AAA
% 4 collagen fiber families, SM, and elastin 
% solid output file   name+"_s%i.dat" %i=0, 1, 2, 3, ....
% for testing code use >> AAA_main(62.01, 0.01, 100, 'test','false',[88.38 1361.1 12.65 0.203])
% Latest update (March 17, 2010) 


function [result] = AAA_main( k_sigma_f, k_sigma_m, Le_low, damage_k,damage_sigma,damage_mu, days, name, parallel)

%global kc P_a r_h H_h nu_e0 nu_f0 nu_m0 phi0 G_h G_e G_m Sa La_M La_0 sigma_f0 sigma_m0 n_elm kq_c kq_m age_max op_time

%addpath('Utility');

if nargin < 5
    parallel = false;
end

%if nargin < 7
damage_params.mu = damage_mu;
damage_params.sigma = damage_sigma;
damage_params.k = damage_k;
%end

format long e;

n_elm=110;             %number of element     30
n_dt = 1/30;             %time steps in a day   1/5
Length = 30;           %longitudinal length =Length *r_h; 
age_max = 300;  % maximum life span of collagen and smooth muscle (days)
kq_c=1/50.0;  % the rate of degradation for collagen
kq_m=1/50.0;  %  .. for smooth muscle

op_time = 210000;

%parameters for constitutive functions (these should be obtained from
%parameter estimation with a healthy aorta
%also note that these parameter with the deposition stretch should satisfy
%the homeostatic stress
%kc=[c1 (elastin), c2(fiber), c3(fiber), c4 (SM), c5 (SM)]
% axis 1 -circumferential direction
% axis 2 -longitudinal direction

rho=1050;  %density of the wall

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters Selection

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Sajjad age 65, SM removed, the chosen one
% kc = [100.63  1686.96  25.865  15.2 11.4]; 
% 
% P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
% r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
% nu_e0 = 0.109;                       %mass fraction of elastin at an initial state
% nu_m0 = 0.15;                      %mass fraction of SM
% nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
%                                    % nu_e0+nu_m0+sum(nu_f0)=1.0
% phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)
% 
% G_h = 1.032 ;                     % homeostatic stretch of fibers
% G_m = 1.1 ;
% G_e = [1.382  1.283];          % initial stretch of elastin layer
% 
% %parameters for active tone
% Sa = 54e+003;  % Max. vasoactive parameter (Pa)
% La_M = 1.4 ;
% La_0 = 0.8 ;

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % initial ones
% kc = [1.4704e+002  3.1804e+003  2.1766e+001  8.3965e+001  5.4094e+000]; 
% 
% P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
% r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
% nu_e0 = 0.2;                       %mass fraction of elastin at an initial state
% nu_m0 = 0.2;                      %mass fraction of SM
% nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
%                                    % nu_e0+nu_m0+sum(nu_f0)=1.0
% phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)
% 
% G_h = 1.04 ;                     % homeostatic stretch of fibers
% G_m = 1.1 ;
% G_e = [1.086  1.0906];          % initial stretch of elastin layer
% 
% %parameters for active tone
% Sa = 36e+003;  % Max. vasoactive parameter (Pa)
% La_M = 1.4 ;
% La_0 = 0.8 ;
% 

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Wilson 2013
% kc = [72  1136  11.2  15.2  11.4]; 
% 
% P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
% r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
% nu_e0 = 0.23;                       %mass fraction of elastin at an initial state
% nu_m0 = 0.15;                      %mass fraction of SM
% nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
%                                    % nu_e0+nu_m0+sum(nu_f0)=1.0
% phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)
% 
% G_h = 1.062 ;                     % homeostatic stretch of fibers
% G_m = 1.1 ;
% G_e = [1.34  1.25];          % initial stretch of elastin layer
% 
% %parameters for active tone
% Sa = 54e+003;  % Max. vasoactive parameter (Pa)
% La_M = 1.4 ;
% La_0 = 0.8 ;

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Shahrokh 2010
% kc = [112  917  25  27  8.5]; 
% 
% P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
% r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
% nu_e0 = 0.23;                       %mass fraction of elastin at an initial state
% nu_m0 = 0.15;                      %mass fraction of SM
% nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
%                                    % nu_e0+nu_m0+sum(nu_f0)=1.0
% phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)
% 
% G_h = 1.07 ;                     % homeostatic stretch of fibers
% G_m = 1.2 ;
% G_e = [1.25  1.25];          % initial stretch of elastin layer
% 
% %parameters for active tone
% Sa = 50e+003;  % Max. vasoactive parameter (Pa)
% La_M = 1.2 ;
% La_0 = 0.7 ;
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Sajjad age 45
% kc = [88.38  1361.12  12.65  73.92  20.15]; 
% 
% P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
% r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
% nu_e0 = 0.203;                       %mass fraction of elastin at an initial state
% nu_m0 = 0.15;                      %mass fraction of SM
% nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
%                                    % nu_e0+nu_m0+sum(nu_f0)=1.0
% phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)
% 
% G_h = 1.088 ;                     % homeostatic stretch of fibers
% G_m = 1.269 ;
% G_e = [1.382  1.283];          % initial stretch of elastin layer
% 
% %parameters for active tone
% Sa = 54e+003;  % Max. vasoactive parameter (Pa)
% La_M = 1.4 ;
% La_0 = 0.8 ;
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Sajjad age 60
% kc = [98.33  1625.87  22.56  73.92  37.10]; 
% 
% P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
% r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
% nu_e0 = 0.132;                       %mass fraction of elastin at an initial state
% nu_m0 = 0.15;                      %mass fraction of SM
% nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
%                                    % nu_e0+nu_m0+sum(nu_f0)=1.0
% phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)
% 
% G_h = 1.042 ;                     % homeostatic stretch of fibers
% G_m = 1.269 ;
% G_e = [1.382  1.283];          % initial stretch of elastin layer
% 
% %parameters for active tone
% Sa = 54e+003;  % Max. vasoactive parameter (Pa)
% La_M = 1.4 ;
% La_0 = 0.8 ;
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Sajjad age 70
% kc = [102.6  1739.33  29.17  73.92  48.4]; 
% 
% P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
% r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
% nu_e0 = 0.085;                       %mass fraction of elastin at an initial state
% nu_m0 = 0.15;                      %mass fraction of SM
% nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
%                                    % nu_e0+nu_m0+sum(nu_f0)=1.0
% phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)
% 
% G_h = 1.023 ;                     % homeostatic stretch of fibers
% G_m = 1.269 ;
% G_e = [1.382  1.283];          % initial stretch of elastin layer
% 
% %parameters for active tone
% Sa = 54e+003;  % Max. vasoactive parameter (Pa)
% La_M = 1.4 ;
% La_0 = 0.8 ;
% % 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Sajjad age 45, SM removed
% kc = [88.38  1361.12  12.65  15.2  11.4]; 
% 
% P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
% r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
% nu_e0 = 0.203;                       %mass fraction of elastin at an initial state
% nu_m0 = 0.15;                      %mass fraction of SM
% nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
%                                    % nu_e0+nu_m0+sum(nu_f0)=1.0
% phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)
% 
% G_h = 1.088 ;                     % homeostatic stretch of fibers
% G_m = 1.1 ;
% G_e = [1.382  1.283];          % initial stretch of elastin layer
% 
% %parameters for active tone
% Sa = 54e+003;  % Max. vasoactive parameter (Pa)
% La_M = 1.4 ;
% La_0 = 0.8 ;
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sajjad age 60, SM removed
kc = [98.33  1625.87  22.56  15.2  11.4]; 

P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
nu_e0 = 0.132;                       %mass fraction of elastin at an initial state
nu_m0 = 0.15;                      %mass fraction of SM
nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
                                   % nu_e0+nu_m0+sum(nu_f0)=1.0
phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)

G_h = 1.042 ;                     % homeostatic stretch of fibers
G_m = 1.1 ;
G_e = [1.382  1.283];          % initial stretch of elastin layer

%parameters for active tone
Sa = 54e+003;  % Max. vasoactive parameter (Pa)
La_M = 1.4 ;
La_0 = 0.8 ;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Sajjad age 70, SM removed
% kc = [102.6  1739.33  29.17  15.2  11.4]; 
% 
% P_a = 100*133.322668;              %average (invivo) pressure : ~100 mmHg       
% r_h= 0.01;                         % invivo mean diameter of the aorta = 2 cm
% nu_e0 = 0.085;                       %mass fraction of elastin at an initial state
% nu_m0 = 0.15;                      %mass fraction of SM
% nu_f0 = [0.1  0.1  0.4  0.4]*(1-nu_e0-nu_m0);  %mass fraction of each fiber at initial state
%                                    % nu_e0+nu_m0+sum(nu_f0)=1.0
% phi0=45.0*pi/180.0;                %alignment angle of diagonal fiber wrt axis 2 (longitudinal)
% 
% G_h = 1.023 ;                     % homeostatic stretch of fibers
% G_m = 1.1 ;
% G_e = [1.382  1.283];          % initial stretch of elastin layer
% 
% %parameters for active tone
% Sa = 54e+003;  % Max. vasoactive parameter (Pa)
% La_M = 1.4 ;
% La_0 = 0.8 ;
% 
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sigma_f0 = 170.35*10^3 ;               % homeostatic stress of collagen (pa)
sigma_m0 = 48.438*10^3 ;               % homeostatic stress for the sum of passive + active


% checking the compatibility conditions between parameters, the
% homeostatic stresses, and balance equation
t_11_c = rho*nu_f0(2)/(1-nu_e0-nu_m0)*dWcdx(G_h, kc)*G_h ...
    + 2* rho*nu_f0(3)/(1-nu_e0-nu_m0)*dWcdx(G_h, kc)*G_h*sin(phi0)^2;
t_22_c = rho*nu_f0(1)/(1-nu_e0-nu_m0)*dWcdx(G_h, kc)*G_h ...
    + 2* rho*nu_f0(3)/(1-nu_e0-nu_m0)*dWcdx(G_h, kc)*G_h*cos(phi0)^2;
t_m = rho*dWmdx(G_m, kc)*G_m+Sa*(1-(La_M-1.0)^2/(La_M-La_0)^2);
t_e =  rho*dWedx(1, G_e(1), G_e(2), kc)*G_e(1);

% Calculation of in vivo thickness

stress = nu_e0*t_e + nu_m0*t_m + (1-nu_e0-nu_m0)*t_11_c;
H_h =  1.0556e-003;

% invivo thickness (m)
%  disp((t_11_c-sigma_f0)/sigma_f0);
%  disp((t_m-sigma_m0)/sigma_m0);
%  disp((H_h-P_a*r_h/stress)/(P_a*r_h/stress));

sigma_f0=t_11_c;
sigma_m0=t_m;
H_h = P_a*r_h/stress;


% Then, initial value at the (invivo) reference configuration are calculated
% at the Gauss points.

%fa_init(name,Data_t0); %without initial damage

if parallel
    % 	% --- Create worker pool
    % 	N = 5;
    % 	poolobj = gcp('nocreate');
    % 	if isempty(poolobj)
    % 		poolsize = 0;
    % 	else
    % 		poolsize = poolobj.NumWorkers;
    % 	end
    %
    % 	if poolsize == 0
    % 		parpool('local',N);
    % 	else
    % 		if poolsize~=N
    % 			delete(poolobj);
    % 			parpool('local',N);
    % 		end
    % 	end
%     
    result = growth_par_elmloop(damage_params, days, k_sigma_f, ...
        k_sigma_m, name, Length, n_dt, floor(100*n_dt), kc, P_a, r_h, H_h, nu_e0, nu_f0, nu_m0, phi0, G_h, G_e, G_m, Sa, La_M, La_0, sigma_f0, sigma_m0, n_elm, kq_c, kq_m, age_max, op_time, Le_low);
    
    % 	delete(poolobj);
else
    % global kc P_a r_h H_h nu_e0 nu_f0 nu_m0 phi0 G_h G_e G_m Sa La_M La_0 sigma_f0 sigma_m0 n_elm kq_c kq_m age_max op_time; %#ok<TLEV>
    growth(days, k_sigma_f, k_sigma_m, name, Length, n_dt, floor(100*n_dt));
end

%clear *

end



