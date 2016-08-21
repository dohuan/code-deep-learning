function [er, estimate] = nntest_regression(nn, x, y)
    estimate = nnpredict(nn, x);
    [~, expected] = max(y,[],2);
    bad = find(estimate ~= expected);    
    er = numel(bad) / size(x, 1);
end
