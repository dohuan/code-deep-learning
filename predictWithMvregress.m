close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
patientID{1} = 'P1'; %G
patientID{2} = 'P2'; %H
patientID{3} = 'P3'; %J
patientID{4} = 'P4'; %K
patientID{5} = 'P5'; %P11
patientID{6} = 'P6'; %P12
patientID{7} = 'P7'; %P13
% --- try mvregress matlab
fileList = dir('./data/real_scans/*.mat');
train = {};
test = {};
traint = [];
testt = [];
for i = 1:length(fileList)
    loadFile = ['./data/real_scans/' fileList(i).name];
    load(loadFile);
    if (length(data)<4)
        fprintf('Discarded: insufficient data\n');
    else
        rawdata = [];
        tmpt = [];
        for j=length(data)-3:length(data)
            if j~=length(data)
                p = polyfit(data(j).pos, data(j).maxd,7);
                x = linspace(min(data(j).pos), max(data(j).pos), 221);
                tmp = [polyval(p,x)];
                rawdata = [rawdata;tmp];
                tmpt = [tmpt, data(j).stime];
            else
                p = polyfit(data(j).pos, data(j).maxd,7);
                x = linspace(min(data(j).pos), max(data(j).pos), 221);
                test_y = [polyval(p,x)];
                testt = [testt;data(j).stime];
            end
        end
        train = [train,rawdata];
        test = [test,test_y];
        traint = [traint;tmpt];
    end
end

nt = size(traint,2);
nd = size(traint,1);
% for i=1:size(traint,1)
%     X{i} = getEstimateMatrix(traint(i,:));
%     x{i} = getEstimateMatrix(testt(i));
% end
est_y = zeros(size(traint,1),221);
for i=1:221
    fprintf('Processing site: %d/221\n',i);
    Y = [];
    for j=1:length(train)
        Y = [Y;(train{j}(:,i))'];
    end
    %model = @(PHI,t)(PHI(:,1))./(1+exp(-(t-PHI(:,2))./PHI(:,3)));
    %model = @(PHI,t)t.*(PHI(:,1)+PHI(:,2).*t+PHI(:,3).*t.^2+PHI(:,4).*t.^3+PHI(:,5).*t.^4);
    model = @(PHI,t)PHI(:,1)+PHI(:,2).*t;
    NUMS = repmat((1:nd)',[1 nt]);
    beta0 = [100 100]; % 5 5 1
    paramorder = length(beta0);
    [beta2,PSI2,stats2,b2] = nlmefit(traint(:),Y(:),...
        NUMS(:),[],model,beta0,'REParamsSelect',[1:paramorder],'RefineBeta0','off');
    %PHI = repmat(beta2,1,nd) + ...          % Fixed effects
    %    [b2(1,:);zeros(1,nd);b2(2,:)];    % Random effects
    PHI = repmat(beta2,1,nd) + ...          % Fixed effects
        [zeros(1,nd); b2(2:end,:)];    % Random effects
    %tmp = repmat(beta2,1,nd);
    %tmp(1,:) = zeros(1,nd);
    %PHI =  tmp + ...          % Fixed effects
        b2;    % Random effects
    for j=1:nd
        %fitted_model=@(t)(PHI(1,j))./(1+exp(-(t-PHI(2,j))./ ...
        %    PHI(3,j)));
        %fitted_model=@(t)t.*(PHI(1,j)+PHI(2,j).*t+PHI(3,j).*t.^2+PHI(4,j).*t.^3+PHI(5,j).*t.^4);
        fitted_model=@(t)PHI(1,j)+PHI(2,j).*t;
        est_y(j,i) = fitted_model(testt(j));
    end
end

err = [];
for i=1:length(test)
    err = [err;rmseCal(est_y(i,:)./max(est_y(i,:)),test{i})];
    fprintf('RMSE (normalized) of patient %s %.2f\n',patientID{i},rmseCal(est_y(i,:)./max(est_y(i,:)),test{i}./max(test{i})));
end                                

fprintf('RMSE: %.2f\n',mean(err));

figure(1)
hold on
for i=1:7
    plot(est_y(i,:));
end
hold off
%save('./mvgress_try')
