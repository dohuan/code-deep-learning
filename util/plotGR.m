function h = plotGR(W)
% W: [hidden/input units]x[visible/output units]
[hid_size,~] = size(W);
s =  sqrt(hid_size);
if (s==floor(s))
    h = figure('name','RBM weights');
    for i=1:hid_size
        subplot(s,s,i)
        plot(W(i,:));
    end
else 
    fprintf('Hidden units are not square\n');
end
end