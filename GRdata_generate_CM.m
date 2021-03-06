close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
addpath(genpath('./Utility'))
%%
tic
load ./data/param_ksigmu

range.k = [0.196 0.287];
range.sigma = [1.72 2.70];
range.mu = [8.77 9.71];

para_info(1).mean = (range.k(1)+range.k(2))/2;
para_info(1).std = (range.k(2)-range.k(1))/2;
para_info(1).name = 'damage_k';
para_info(1).dist_type = 'Normal';
para_info(1).trunc = range.k;

para_info(2).mean = (range.sigma(1)+range.sigma(2))/2;
para_info(2).std = (range.sigma(2)-range.sigma(1))/2;
para_info(2).name = 'damage_sigma';
para_info(2).dist_type = 'Normal';
para_info(2).trunc = range.sigma;

para_info(3).mean = (range.mu(1)+range.mu(2))/2;
para_info(3).std = (range.mu(2)-range.mu(1))/2;
para_info(3).name = 'damage_mu';
para_info(3).dist_type = 'Normal';
para_info(3).trunc = range.mu;

[coeff, Ortho] = getCM(para_info,1);
poly_order = length(Ortho);

y = [];
for i=1:size(param,1)
    H_tmp = [];
    for k1=1:poly_order+1
        for k2=1:poly_order+1
            for k3=1:poly_order+1
                H_tmp = [H_tmp,Ortho{1}.H{k1}(param(i,1))*Ortho{2}.H{k2}(param(i,2))*...
                    Ortho{3}.H{k3}(param(i,3))];
            end
        end
    end
    y = [y; H_tmp*coeff];
end

% --- chop data into feature and target then save
data.train_x = [];
data.train_y = [];
fnum = 3; % take [fnum] scans in the past as features, predict the scan [fnum+1]-th
t_offset = 12*(fnum+1); % t_offset = 12(months)*(fnum+1)(scans)*30(days)/30(days)
t_step = 12;
for i=1:size(y,1)
    rawdata = (reshape(y(i,:),221,[]))';
    ran_ix = randi([1,size(rawdata,1)-t_offset],100,1);
    for j=1:length(ran_ix)
        tmp_x = [];
        stoken = (rand()>.7); % shift token
        for k=1:fnum+1
            ix = ran_ix(j)+t_step*(k-1);
            if (stoken ==1 )
                data_tmp = (curve_shifter(single(rawdata(ix,:)'), .7))';
            else
                data_tmp = single(rawdata(ix,:));
            end
            if k~=fnum+1
                tmp_x = [tmp_x, data_tmp];
            else
                
                tmp_y = data_tmp;
            end
        end
        data.train_x = [data.train_x;tmp_x];
        data.train_y = [data.train_y; tmp_y];
    end
end

% --- create data with 2-peak AAA
for i=1:4:size(y,1)
    rawdata = (reshape(y(i,:),221,[]))';
    ran_ix1 = randi([1,size(rawdata,1)-t_offset],1,1);
    ran_ix2 = randi([1,size(rawdata,1)-t_offset],1,1);
    for j=1:length(ran_ix1)
        tmp_x = [];
        %stoken = (rand()>.7); % shift token
        for k=1:fnum+1
            ix1 = ran_ix1(j)+t_step*(k-1);
            ix2 = ran_ix2(j)+t_step*(k-1);
            ratio = .05+(.1-.05).*rand();
            
            data_tmp = curve_merger(rawdata(ix1,:),rawdata(ix2,:),1,ratio);
            data_tmp = single(data_tmp);
            if k~=fnum+1
                tmp_x = [tmp_x, data_tmp];
            else
                
                tmp_y = data_tmp;
            end
        end
        data.train_x = [data.train_x;tmp_x];
        data.train_y = [data.train_y; tmp_y];
    end
    
end

% --- mix all data together
mixix = randperm(size(data.train_x,1));
data.train_x = data.train_x(mixix,:);
data.train_y = data.train_y(mixix,:);

data.fnum = fnum;
save('./data/dataCM-augmented','data');
time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);






