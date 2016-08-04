function dbn = dbntrain(dbn, x, opts)
n = numel(dbn.rbm);

fprintf('Traning RBM: 1\n');
if (strcmp(opts.visibleDist,'Gauss')==1)
    dbn.rbm{1} = rbmtrainGauss(dbn.rbm{1}, x, opts);
else
    dbn.rbm{1} = rbmtrain(dbn.rbm{1}, x, opts);
end

for i = 2 : n
    fprintf('Traning RBM layer: %d\n',i);
    x = rbmup(dbn.rbm{i - 1}, x);
    dbn.rbm{i} = rbmtrain(dbn.rbm{i}, x, opts);
end

end
