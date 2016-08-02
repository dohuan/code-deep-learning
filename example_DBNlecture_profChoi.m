close all
clear
clc
set(0,'defaultfigurecolor',[1 1 1])
%%
tic
load ./data/mnist_uint8;

train_x = double(train_x) / 255;
test_x  = double(test_x)  / 255;
train_y = double(train_y);
test_y  = double(test_y);

%%  ex1 train a 100 hidden unit RBM and visualize its weights

% --- 1 layer of hidden unit with size 100
%dbn.sizes = 100;
% --- 2 layers of hidden unit with size 100
dbn.sizes = [100 100];

opts.numepochs =   1;
opts.batchsize = 100;
opts.momentum  =   0;
opts.alpha     =   1;
dbn = dbnsetup(dbn, train_x, opts);
dbn = dbntrain(dbn, train_x, opts);
figure(1)
for i=1:length(dbn.sizes)-1
    subplot(1,length(dbn.sizes)-1,i)
    a = visualize(dbn.rbm{i}.W');   %  Visualize the RBM weights
    imagesc(a, [min(a(:)) max(a(:))]);
    axis tight
    colormap gray
    if i==length(dbn.sizes)-1
        title(sprintf('RBM layer: %d (most general)',i))
    else
        title(sprintf('RBM layer: %d',i))
    end
end

%% Use pre-trained weights from DBN to form a Neural Network
% --- unfold dbn to nn
nn = dbnunfoldtonn(dbn, 10);
nn.activation_function = 'sigm';

% --- train nn
opts.numepochs =  3;
opts.batchsize = 100;
nn = nntrain(nn, train_x, train_y, opts);
[er, bad, est_label] = nntest(nn, test_x, test_y);

% --- Plot 10 wrong cases
figure(2)
for i =1:10
    badIx = bad(i);
    label = find(test_y(badIx,:)==1)-1;
    h = rot90(reshape(test_x(badIx,:),28,28));
    subplot(2,5,i)
    pcolor(h)
    axis tight
    box on
    title(sprintf('Est: %d Label: %d',est_label(badIx)-1, label))
end


time_run = toc;
fprintf('\nRun time: %.2f minutes\n',time_run/60);

