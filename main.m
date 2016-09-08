close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
tStart = tic;
patientID{1} = 'G';
patientID{2} = 'H';
patientID{3} = 'J';
patientID{4} = 'K';
patientID{5} = 'P11';
patientID{6} = 'P12';
patientID{7} = 'P13';

% load ./data/dataGR-augmented;
% fnum = 3;
% for i=1:size(data.train_x,1)
%     train_x(i,:) = data.train_x(i,:)./max(data.train_x(i,:));
% end
% for i=1:size(data.train_y,1)
%     train_y(i,:) = data.train_y(i,:)./max(data.train_y(i,:));
%     %train_y(i,:) = data.train_y(i,:);
% end
% clear data

%load ./data/data_train
load ./data/dataGR-augmented
train_x = double(data.train_x);
train_y = double(data.train_y);
clear data


%%  ex1 train a 100 hidden unit RBM and visualize its weights

% --- 1 layer of hidden unit with size 100
%dbn.sizes = 100;
% --- 2 layers of hidden unit with size 100
dbn.sizes = [100 100];

opts.numepochs =   3;
opts.batchsize = 100;
opts.momentum  =   0;
opts.alpha     =   5E-3;   % alpha: learn rate
opts.visibleDist   = 'Gauss'; % 'Gauss' or 'binomial'
dbn = dbnsetup(dbn, train_x, opts);
dbn = dbntrain(dbn, train_x, opts);

% --- Plot weights (comment if not necesary)
% figure(1)
% for i=1:length(dbn.sizes)-1
%     plotGR(dbn.rbm{i}.W);
%     if i==length(dbn.sizes)-1
%         title(sprintf('RBM layer: %d (most general)',i))
%     else
%         title(sprintf('RBM layer: %d',i))
%     end
% end

%% Use pre-trained weights from DBN to form a Neural Network
% --- unfold dbn to nn
%nn = dbnunfoldtonn(dbn, 10);
nn = dbnunfoldtonn(dbn, size(train_y,2));
nn.activation_function = 'sigm';

% --- train NN USING GR data
opts.numepochs =  1;
opts.batchsize = 100;
nn.learningRate = .03;
nn = nntrain(nn, train_x, train_y, opts);

% --- train NN USING REAL data (for better fine tuning)
%load ./data/dataREAL_normalized_uni_max
load ./data/dataREAL_unnormalized
opts.numepochs =  1;
opts.batchsize = 1;
nn.learningRate = .0003;
nn = nntrain(nn, data.ft_x, data.ft_y, opts);

est = nnpredict(nn,data.test_x);

figure(2)
hold on
plot(est(1,:))
plot(est(2,:))

figure(1)
for i=1:size(data.test_y,1)
    subplot(2,4,i)
    hold on
    plot(est(i,:),'LineWidth',2)
    plot(data.test_y(i,:),'LineWidth',2)
    hold off
    legend('estimated','true')
    title(patientID{i})
end





time_run = toc(tStart);
fprintf('\nRun time: %.2f minutes\n',time_run/60);