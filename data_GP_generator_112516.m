close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
addpath(genpath('./code-deep-learning/gpml'));
fileList = dir('./data/real_scans/*.mat');

X = [];
covfunc = @covSEard;
likfunc = @likGauss;
meanfunc = @meanZero;
fitOrder = 5;

for i =1:length(fileList)
    loadFile = ['./data/real_scans/' fileList(i).name];
	load(loadFile);
    max_dia = [];
    time_stamp = [];
    for j=1:length(data)
        p = polyfit(data(j).pos, data(j).maxd,fitOrder);
        x = linspace(min(data(j).pos), max(data(j).pos), 221);
        tmp = polyval(p,x);
        max_dia=[max_dia;tmp];
        time_stamp = [time_stamp;data(j).time];
    end
    meanMaxDia = mean(max_dia,1);
    for j=1:size(max_dia,1)
        max_dia(j,:)=max_dia(j,:)-meanMaxDia;
    end
    timerange = (round(min(time_stamp)):1:round(max(time_stamp)))';
    
    Xtmp = zeros(length(timerange),size(max_dia,2));
    
    
    % ---
    for j=1:size(max_dia,2)
        hyp.cov(1) = log(3.2);
        hyp.cov(2) = log(20);
        hyp.lik = log(0.03);
        hyp = minimize(hyp, @gp, -100, @infExact, [], ...
            covfunc, likfunc, time_stamp,max_dia(:,j));
        
        [Xtmp(:,j), ~] = gp(hyp, @infExact, [], covfunc, likfunc, ...
                time_stamp,max_dia(:,j), timerange);
        Xtmp(:,j) = Xtmp(:,j) + meanMaxDia(j)*ones(size(Xtmp,1),1);
        
    end
    
    % ---
    
    X = [X;Xtmp];
end


hold on
for i=1:size(Xtmp,1)
plot(Xtmp(i,:));
end
axis tight