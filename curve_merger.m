function out = curve_merger(curve1,curve2,ratio1,ratio2)
if (length(curve1)~=length(curve2))
    error('Lengths of the two curves are not the same.\n');
end
curve_tmp = fliplr((curve_shifter(curve2', .6))');
offset = min(curve1);
curve_tmp = curve_tmp - offset;
%ratio_ = ratio2*max(curve2)/max(curve1);
out = curve1.*ratio1+curve_tmp.*ratio2;
end