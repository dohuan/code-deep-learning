close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
addpath(genpath('./code-deep-learning/gpml'));
fileList = dir('./data/real_scans/*.mat');

covfunc = @covSEard;
likfunc = @likGauss;
meanfunc = @meanZero;

hyper = [];
numRan = 10;
numObs = 100;
for i =1:length(fileList)
    loadFile = ['./data/real_scans/' fileList(i).name];
    load(loadFile);
    for j=1:length(data)
        feature = data(j).pos;
        feature = feature-mean(feature);
        target = data(j).maxd;
        nt = length(target);
        
        for k=1:numRan
            ix = randperm(nt);
            ix = ix(1:numObs);
            X = feature(ix,:);
            y = target(ix);
            
            hyp.cov(1) = log(1);
            hyp.cov(2) = log(1);
            hyp.lik = log(0.03);
            hyp = minimize(hyp, @gp, -500, @infExact, [], ...
                covfunc, likfunc, X,y);
            
            hyper = [hyper;[exp(hyp.cov(1)) exp(hyp.cov(2))]];
        end
    end
end

figure(1)
subplot(1,2,1)
hist(hyper(:,1),100);
title('scale')
subplot(1,2,2)
hist(hyper(:,2),100);
title('feature bandwidth')