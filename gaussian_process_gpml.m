function [mu, sig] = gaussian_process_gpml(x,y,x_test)
sampleRatio = 100;
sampleIx = randi([1, size(x,1)],round(size(x,1)/sampleRatio),1);
%sampleIx = 1:sampleRatio:size(x,1);

x = x(sampleIx,:);
y = y(sampleIx,:);

covfunc = @covSEard;
likfunc = @likGauss;

hyp.cov(1) = log(20.2); % .2
hyp.cov(2) = log(3.2); % .2
hyp.cov(3) = log(3.2); % .2
hyp.cov(4) = log(3.2); % .2
hyp.lik = log(0.05);   % .05

hyp = minimize(hyp, @gp, -30, @infExact, [], covfunc, likfunc, x, y);
[mu, sig] = gp(hyp, @infExact, [], covfunc, likfunc, x, y, x_test);
end