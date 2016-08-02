function dbn = dbntrain(dbn, x, opts)
    n = numel(dbn.rbm);
    
    fprintf('Traning RBM: 1\n');
    dbn.rbm{1} = rbmtrain(dbn.rbm{1}, x, opts);
    for i = 2 : n
        fprintf('Traning RBM: %d\n',i);
        x = rbmup(dbn.rbm{i - 1}, x);
        dbn.rbm{i} = rbmtrain(dbn.rbm{i}, x, opts);
    end

end
