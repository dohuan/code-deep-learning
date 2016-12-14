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
    feature = [];
    time_stamp = [];
    xmin = 1000;
    xmax = 0;
    for j=1:length(data)
        p = polyfit(data(j).pos, data(j).maxd,fitOrder);
        %x = linspace(min(data(j).pos), max(data(j).pos), 221);
        %tmp = polyval(p,x);
        %max_dia=[max_dia;tmp];
        if (max(data(j).pos)>xmax)
            xmax = max(data(j).pos);
        end
        if (min(data(j).pos)<xmin)
            xmin = min(data(j).pos);
        end
        feature = [feature;p];
        time_stamp = [time_stamp;data(j).time];
    end
    
    xrange = linspace(xmin,xmax,221);
    meanfeature = mean(feature,1);
    feature_ = feature;
    for j=1:size(feature,1)
        feature(j,:)=feature(j,:)-meanfeature;
    end
    timerange = (round(min(time_stamp)):1:round(max(time_stamp)))';
    
    Xtmp = zeros(length(timerange),length(xrange));
    
    
    % ---
    p_ = [];
    for j=1:size(feature,2)
        hyp.cov(1) = log(5.2);
        hyp.cov(2) = log(15.5);
        hyp.lik = log(0.03);
        hyp = minimize(hyp, @gp, -100, @infExact, [], ...
            covfunc, likfunc, time_stamp,feature(:,j));
        
        [ptmp, ~] = gp(hyp, @infExact, [], covfunc, likfunc, ...
                time_stamp,feature(:,j), timerange);
        ptmp = ptmp + meanfeature(j)*ones(size(Xtmp,1),1);
        p_ = [p_, ptmp];
    end
    
    for j=1:size(Xtmp,1)
        Xtmp(j,:) = polyval(p_(j,:),xrange);
    end
    
    % ---
    
    X = [X;Xtmp];
end


hold on
for i=1:size(Xtmp,1)
plot(Xtmp(i,:));
end
for i=1:length(data)
    p = polyfit(data(i).pos, data(i).maxd,fitOrder);
    plot(polyval(p,xrange),'ko-');
end
axis tight



hold on
for i=1:size(p_,1)
plot(p_(i,:));
end
for i=1:size(feature_,1)
plot(feature_(i,:),'ko');
end
axis tight