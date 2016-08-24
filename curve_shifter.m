function out = curve_shifter(x, ratio)
%% shift a curve in a Mobius manner
%--- find axial location where ratio happens
%--- x: nx1 vector
endvalue = x(end);
cutvalue = (max(x(:))-min(x(:)))*ratio + min(x(:));
distance = (x-cutvalue).^2;
[~,ix] = sort(distance,'ascend');
cutix = min(ix(1:2));
out = x(cutix:end);
gap = length(x)-length(out);
if (gap<=0)
    fprintf('Something is wrong!\n');
else
    out = [out; ones(gap,1)*endvalue];
end

end