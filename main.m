close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
tic
load ./data/dataGR-full;

for i=1:size(data.train_x,1)
    train_x(i,:) = data.train_x(i,:)./max(data.train_x(i,:));
end
for i=1:size(data.train_y,1)
    train_y(i,:) = data.train_y(i,:)./max(data.train_y(i,:));
    %train_y(i,:) = data.train_y(i,:);
end
clear data

% --- Load real data
load ./data/BC_20160608/H.mat
test_x = [];
test_y = [];
for i=length(H)-4:1:length(H)
    p = polyfit(H(i).pos, H(i).maxd,7);
    x = linspace(min(H(i).pos), max(H(i).pos), 221);
    tmp = polyval(p,x);
    if (i~=length(H))
        test_x = [test_x, tmp];
    else
        test_y = tmp./max(tmp);
    end
end
test_x = test_x./(max(test_x));



%%  ex1 train a 100 hidden unit RBM and visualize its weights

% --- 1 layer of hidden unit with size 100
%dbn.sizes = 100;
% --- 2 layers of hidden unit with size 100
dbn.sizes = [100 100 100];

opts.numepochs =   1;
opts.batchsize = 100;
opts.momentum  =   0;
opts.alpha     =   1;
opts.visibleDist   = 'bino'; % 'Gauss' or 'binomial'
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

% --- train nn
opts.numepochs =  1;
opts.batchsize = 100;
nn = nntrain(nn, train_x, train_y, opts);
[er, bad, est_label] = nntest(nn, test_x, test_y);

figure(1)
hold on
plot(est_label,'LineWidth',2)
plot(test_y,'LineWidth',2)
hold off
legend('estimated','true')

time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);