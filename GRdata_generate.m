close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
addpath(genpath('./Utility'))
%%
%range.k = [0.196 0.287];
%range.sigma = [1.72 2.70];
%range.mu = [8.77 9.71];
tic
load param_ksigmu
k_sigma_f = 0.05;
k_sigma_m = 0.01;
Le_lowH = 12.027458;             % patient H
Le_lowK = 1.302589479770089e+01; % patient K


runPC = 'ME197'; % ME197 or DECS

if (strcmp(runPC,'ME197')==1)
    run_param = param(1:256,:);
    for i=1:size(run_param,1)
        fprintf('Running on %s iteration: %d\n',runPC,i);
        damage_k        = run_param(i,1);
        damage_sigma    = run_param(i,2);
        damage_mu       = run_param(i,3);
        result = AAA_main(k_sigma_f,k_sigma_m,Le_lowH, damage_k,...
                                  damage_sigma,damage_mu, 5000,'H','true');
        savefile = ['GR_H_' runPC '_' num2str(i)];
        save(savefile);
    end
else
    run_param = param(257:end,:);
    for i=1:size(run_param,1)
        fprintf('Running on %s iteration: %d\n',runPC,i);
        damage_k        = run_param(i,1);
        damage_sigma    = run_param(i,2);
        damage_mu       = run_param(i,3);
        result = AAA_main(k_sigma_f,k_sigma_m,Le_lowH, damage_k,...
                                  damage_sigma,damage_mu, 5000,'H','true');
        savefile = ['GR_H_' runPC '_' num2str(i)];
        save(savefile);
    end
end
%damage_k = 0.2363978284020520;
%damage_sigma = 1.819552038180287;
%damage_mu = 12.117;
%savefile = ['GR_H_' runPC];
%save(savefile);


time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);