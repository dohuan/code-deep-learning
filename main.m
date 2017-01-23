close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
addpath(genpath('./code-deep-learning/gpml'));
tStart = tic;
patientID{1} = 'P1'; %G
patientID{2} = 'P2'; %H
patientID{3} = 'P3'; %J
patientID{4} = 'P4'; %K
patientID{5} = 'P5'; %P11
patientID{6} = 'P6'; %P12
patientID{7} = 'P7'; %P13


fnum = 3;
ifLoad = 1; % Run only ONCE then load saved data
if ifLoad==0
    %load ./data/dataGR-augmented
    load ./data/dataCM-augmented
    scale_feature = [];
    for i=1:size(data.train_x,1)
        [ic,~] = sort(data.train_x(i,:),'descend');
        scale_feature = [scale_feature; ic(fnum:-1:1)]; 
        train_x(i,:) = data.train_x(i,:)./ic(1);
    end
    scale_label = [];
    for i=1:size(data.train_y,1)
        train_y(i,:) = data.train_y(i,:)./max(data.train_y(i,:));
        scale_label = [scale_label; max(data.train_y(i,:))];
    end
    clear data
    
    load ./data/dataREAL_unnormalized
    for i=1:size(data.ft_x,1)
        [ic,~] = sort(data.ft_x(i,:),'descend');
        scale_feature = [scale_feature; ic(fnum:-1:1)]; 
        ft_x(i,:) = data.ft_x(i,:)./ic(1);
    end
    for i=1:size(data.ft_y,1)
        ft_y(i,:) = data.ft_y(i,:)./max(data.ft_y(i,:));
        scale_label = [scale_label; max(data.ft_y(i,:))];
    end
    % --- save scale in order: train-finetune
    scaleTrackTrain = [scale_feature scale_label];
    
    
    scale_feature = [];
    scale_label = [];
    for i=1:size(data.test_x,1)
        [ic,~] = sort(data.test_x(i,:),'descend');
        scale_feature = [scale_feature; ic(fnum:-1:1)]; 
        test_x(i,:) = data.test_x(i,:)./ic(1);
    end
    for i=1:size(data.test_y,1)
        test_y(i,:) = data.test_y(i,:)./max(data.test_y(i,:));
        scale_label = [scale_label; max(data.test_y(i,:))];
    end
    
    scaleTrackTest = [scale_feature scale_label];
    
    save('./data/data_all_GR','train_x','train_y','scaleTrackTrain',...
        'scaleTrackTest','ft_x','ft_y','test_x','test_y');
else
    load ./data/data_all_GR
end



%%  ex1 train a 100 hidden unit RBM and visualize its weights
DBNtime = tic;
% --- 1 layer of hidden unit with size 100
%dbn.sizes = 100;
% --- 2 layers of hidden unit with size 100
dbn.sizes = [400 100]; % 100 100 | 900 36 | 400 100(published)

opts.numepochs =   3;
opts.batchsize = 100;
opts.momentum  =   0;
opts.alpha     =   50E-4;   % alpha: learn rate, use 1E-4 for GR data 
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
nn.learningRate = .1;
nn = nntrain(nn, train_x, train_y, opts);

% --- train NN USING REAL data
opts.numepochs =  260;  % 250 (published)
opts.batchsize = 1;
nn.learningRate = .1;
nn = nntrain(nn, ft_x, ft_y, opts);

DBNtimerun = toc(DBNtime);
fprintf('\nDBN training time: %.2f seconds\n',DBNtimerun);

est = nnpredict(nn,test_x);

% --- Use GP here to predict scale for test
[scale_test, ~] = gaussian_process_gpml(scaleTrackTrain(:,1:3),scaleTrackTrain(:,end),...
                                            scaleTrackTest(:,1:3));
err = [];
for i=1:size(test_y,1)
    estUnscaled = est(i,:).*scale_test(i);
    err = [err;rmseCal(estUnscaled,test_y(i,:)*scaleTrackTest(i,end))];
end                                

fprintf('RMSE: %.2f\n',mean(err));

figure(1)
title('unscaled')
for i=1:size(test_y,1)
    subplot(2,4,i)
    hold on
    plot(est(i,:),'LineWidth',2)
    plot(test_y(i,:),'LineWidth',2)
    hold off
    legend('estimated','true')
    title(patientID{i})
end

figure(2)
for i=1:size(test_y,1)
    subplot(2,4,i)
    hold on
    estPlot = est(i,:).*scale_test(i);
    estPlot = smooth(estPlot,.1,'lowess');
    plot(est(i,:).*scale_test(i),'g-.','LineWidth',2)
    plot(test_y(i,:)*scaleTrackTest(i,end),'r:','LineWidth',2)
    plot(estPlot,'b-','LineWidth',2);
    hold off
    legend('prediction','true','smoothed prediction')
    title(patientID{i})
    box on
    axis tight
end

figure(4)
hold on
for i=1:size(est,1)
    plot(est(i,:));
end




time_run = toc(tStart);
fprintf('\nRun time: %.2f minutes\n',time_run/60);