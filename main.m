close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
tic
load ./data/dataGR;

for i=1:size(data.train,1)
    train_x(i,:) = data.train(i,:)./max(data.train(i,:));
end
for i=1:size(data.test,1)
    test_x(i,:) = data.test(i,:)./max(data.test(i,:));
end
clear data

%%  ex1 train a 100 hidden unit RBM and visualize its weights

% --- 1 layer of hidden unit with size 100
%dbn.sizes = 100;
% --- 2 layers of hidden unit with size 100
dbn.sizes = [400 400];

opts.numepochs =   1;
opts.batchsize = 50;
opts.momentum  =   0;
opts.alpha     =   1;
opts.visibleDist   = 'bino'; % 'Gauss' or 'binomial'
dbn = dbnsetup(dbn, train_x, opts);
dbn = dbntrain(dbn, train_x, opts);
figure(1)

for i=1:length(dbn.sizes)-1
    plotGR(dbn.rbm{i}.W);
    if i==length(dbn.sizes)-1
        title(sprintf('RBM layer: %d (most general)',i))
    else
        title(sprintf('RBM layer: %d',i))
    end
end

time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);