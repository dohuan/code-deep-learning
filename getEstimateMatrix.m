function M=getEstimateMatrix(t)
M=zeros(length(t),3);
for i=1:length(t)
    M(i,:) = [1 t(i) t(i).^2];
end
end