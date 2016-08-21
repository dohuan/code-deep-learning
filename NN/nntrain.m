function [nn, L]  = nntrain(nn, train_x, train_y, opts, val_x, val_y)
%NNTRAIN trains a neural net
% [nn, L] = nnff(nn, x, y, opts) trains the neural network nn with input x and
% output y for opts.numepochs epochs, with minibatches of size
% opts.batchsize. Returns a neural network nn with updated activations,
% errors, weights and biases, (nn.a, nn.e, nn.W, nn.b) and L, the sum
% squared error for each training minibatch.

assert(isfloat(train_x), 'train_x must be a float');
assert(nargin == 4 || nargin == 6,'number ofinput arguments must be 4 or 6')

loss.train.e               = [];
loss.train.e_frac          = [];
loss.val.e                 = [];
loss.val.e_frac            = [];
opts.validation = 0;
if nargin == 6
    opts.validation = 1;
end

fhandle = [];
if isfield(opts,'plot') && opts.plot == 1
    fhandle = figure();
end

m = size(train_x, 1);

batchsize = opts.batchsize;
numepochs = opts.numepochs;

numbatches = m / batchsize;

assert(rem(numbatches, 1) == 0, 'numbatches must be a integer');

L = zeros(numepochs*numbatches,1);
n = 1;
for i = 1 : numepochs
    tic;
    
    kk = randperm(m);
    for l = 1 : numbatches
        batch_x = train_x(kk((l - 1) * batchsize + 1 : l * batchsize), :);
        
        %Add noise to input (for use in denoising autoencoder)
        if(nn.inputZeroMaskedFraction ~= 0)
            batch_x = batch_x.*(rand(size(batch_x))>nn.inputZeroMaskedFraction);
        end
        
        batch_y = train_y(kk((l - 1) * batchsize + 1 : l * batchsize), :);
        
        nn = nnff(nn, batch_x, batch_y);
        nn = nnbp(nn);
        nn = nnapplygrads(nn);
        
        L(n) = nn.L;
        
        n = n + 1;
    end
    
    t = toc;

    if opts.validation == 1
        loss = nneval(nn, loss, train_x, train_y, val_x, val_y);
        str_perf = sprintf('; Full-batch train mse = %f, val mse = %f', loss.train.e(end), loss.val.e(end));
    else
        loss = nneval(nn, loss, train_x, train_y);
        str_perf = sprintf('; Full-batch train err = %f', loss.train.e(end));
    end
    if ishandle(fhandle)
        nnupdatefigures(nn, fhandle, loss, opts, i);
    end
        
    disp(['epoch ' num2str(i) '/' num2str(opts.numepochs) '. Took ' num2str(t) ' seconds' '. Mini-batch mean squared error on training set is ' num2str(mean(L((n-numbatches):(n-1)))) str_perf]);
    nn.learningRate = nn.learningRate * nn.scaling_learningRate;
end
end


% function nn = nnff(nn, x, y)
% %NNFF performs a feedforward pass
% % nn = nnff(nn, x, y) returns an neural network structure with updated
% % layer activations, error and loss (nn.a, nn.e and nn.L)
% 
%     n = nn.n;
%     m = size(x, 1);
%     
%     x = [ones(m,1) x];
%     nn.a{1} = x;
% 
%     %feedforward pass
%     for i = 2 : n-1
%         switch nn.activation_function 
%             case 'sigm'
%                 % Calculate the unit's outputs (including the bias term)
%                 nn.a{i} = sigm(nn.a{i - 1} * nn.W{i - 1}');
%             case 'tanh_opt'
%                 nn.a{i} = tanh_opt(nn.a{i - 1} * nn.W{i - 1}');
%         end
%         
%         %dropout
%         if(nn.dropoutFraction > 0)
%             if(nn.testing)
%                 nn.a{i} = nn.a{i}.*(1 - nn.dropoutFraction);
%             else
%                 nn.dropOutMask{i} = (rand(size(nn.a{i}))>nn.dropoutFraction);
%                 nn.a{i} = nn.a{i}.*nn.dropOutMask{i};
%             end
%         end
%         
%         %calculate running exponential activations for use with sparsity
%         if(nn.nonSparsityPenalty>0)
%             nn.p{i} = 0.99 * nn.p{i} + 0.01 * mean(nn.a{i}, 1);
%         end
%         
%         %Add the bias term
%         nn.a{i} = [ones(m,1) nn.a{i}];
%     end
%     switch nn.output 
%         case 'sigm'
%             nn.a{n} = sigm(nn.a{n - 1} * nn.W{n - 1}');
%         case 'linear'
%             nn.a{n} = nn.a{n - 1} * nn.W{n - 1}';
%         case 'softmax'
%             nn.a{n} = nn.a{n - 1} * nn.W{n - 1}';
%             nn.a{n} = exp(bsxfun(@minus, nn.a{n}, max(nn.a{n},[],2)));
%             nn.a{n} = bsxfun(@rdivide, nn.a{n}, sum(nn.a{n}, 2)); 
%     end
% 
%     %error and loss
%     nn.e = y - nn.a{n};
%     
%     switch nn.output
%         case {'sigm', 'linear'}
%             nn.L = 1/2 * sum(sum(nn.e .^ 2)) / m; 
%         case 'softmax'
%             nn.L = -sum(sum(y .* log(nn.a{n}))) / m;
%     end
% end


function nn = nnbp(nn)
%NNBP performs backpropagation
% nn = nnbp(nn) returns an neural network structure with updated weights 
    
    n = nn.n;
    sparsityError = 0;
    switch nn.output
        case 'sigm'
            d{n} = - nn.e .* (nn.a{n} .* (1 - nn.a{n})); 
            % NOTE: here nn.a{n} is the same as h(a) in Bishop 5.58
            % for sigmod: h(a)' = h(a)*(1-h(a))
        case {'softmax','linear'}
            d{n} = - nn.e;
    end
    for i = (n - 1) : -1 : 2
        % Derivative of the activation function
        switch nn.activation_function 
            case 'sigm'
                d_act = nn.a{i} .* (1 - nn.a{i});
            case 'tanh_opt'
                d_act = 1.7159 * 2/3 * (1 - 1/(1.7159)^2 * nn.a{i}.^2);
        end
        
        if(nn.nonSparsityPenalty>0)
            pi = repmat(nn.p{i}, size(nn.a{i}, 1), 1);
            sparsityError = [zeros(size(nn.a{i},1),1) nn.nonSparsityPenalty * (-nn.sparsityTarget ./ pi + (1 - nn.sparsityTarget) ./ (1 - pi))];
        end
        
        % Backpropagate first derivatives
        if i+1==n % in this case in d{n} there is not the bias term to be removed             
            d{i} = (d{i + 1} * nn.W{i} + sparsityError) .* d_act; % Bishop (5.56)
        else % in this case in d{i} the bias term has to be removed
            d{i} = (d{i + 1}(:,2:end) * nn.W{i} + sparsityError) .* d_act;
        end
        
        if(nn.dropoutFraction>0)
            d{i} = d{i} .* [ones(size(d{i},1),1) nn.dropOutMask{i}];
        end

    end

    for i = 1 : (n - 1)
        if i+1==n
            nn.dW{i} = (d{i + 1}' * nn.a{i}) / size(d{i + 1}, 1);
        else
            nn.dW{i} = (d{i + 1}(:,2:end)' * nn.a{i}) / size(d{i + 1}, 1);      
        end
    end
end


function nn = nnapplygrads(nn)
%NNAPPLYGRADS updates weights and biases with calculated gradients
% nn = nnapplygrads(nn) returns an neural network structure with updated
% weights and biases
    
    for i = 1 : (nn.n - 1)
        if(nn.weightPenaltyL2>0)
            dW = nn.dW{i} + nn.weightPenaltyL2 * [zeros(size(nn.W{i},1),1) nn.W{i}(:,2:end)];
        else
            dW = nn.dW{i};
        end
        
        dW = nn.learningRate * dW;
        
        if(nn.momentum>0)
            nn.vW{i} = nn.momentum*nn.vW{i} + dW;
            dW = nn.vW{i};
        end
            
        nn.W{i} = nn.W{i} - dW;
    end
end


