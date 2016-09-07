function x = rbmup(rbm, x)
    x = sigm(repmat(rbm.c', size(x, 1), 1) + x * rbm.W');
%     tmp = (repmat(rbm.c', size(x, 1), 1) + x * rbm.W');
%     for i=1:size(tmp,1)
%         tmp(i,:) = tmp(i,:)-mean(tmp(i,:));
%     end
%     x = sigm(tmp);
end
