function [out,root_poly] = ortho_poly(order,dinfo)
%%     Return orthogonal polynomial for arbitrary normal distribution
% dist_type.name dist_type.mean dist_type.std
syms x

n = order+1;

H{1} = @(x)x.*0;
h{1} = @(x)x.*0;

H{2} = @(x)x.*0+1;
h{2} = @(x)x.*0+1;

for i=3:n+2
    if(i==n+2)
        temp = matlabFunction(x.*h{i-1}(x));
        root_poly.H = matlabFunction(x.*h{i-1}(x)-ortho_int(@(x)temp(x),@(x)h{i-1}(x),dinfo)*h{i-1}(x)-...
            sqrt(ortho_int(@(x)H{i-1}(x),@(x)H{i-1}(x),dinfo))*h{i-2}(x));
        root_poly.h = matlabFunction(root_poly.H(x)/sqrt(ortho_int(@(x)root_poly.H(x),@(x)root_poly.H(x),dinfo)));
    else
        temp = matlabFunction(x.*h{i-1}(x));
        H{i} = matlabFunction(x.*h{i-1}(x)-ortho_int(@(x)temp(x),@(x)h{i-1}(x),dinfo)*h{i-1}(x)-...
            sqrt(ortho_int(@(x)H{i-1}(x),@(x)H{i-1}(x),dinfo))*h{i-2}(x));
        h{i} = matlabFunction(H{i}(x)/sqrt(ortho_int(@(x)H{i}(x),@(x)H{i}(x),dinfo)));
    end
end

% for i=2:order+2
%     if (i==1)
%         H{i} = matlabFunction(x-ortho_int(@(x)x,@(x)x.*0+1,dinfo));
%         h{i} = matlabFunction(H{i}(x)/sqrt(ortho_int(@(x)H{i}(x),@(x)H{i}(x),dinfo)));
%     elseif(i==2)
%         temp = matlabFunction(x.*h{1}(x));
%         H{i} = matlabFunction(x.*h{1}(x)-ortho_int(@(x)temp(x),@(x)h{1}(x),dinfo)*h{1}(x)-...
%             sqrt(ortho_int(@(x)H{1}(x),@(x)H{1}(x),dinfo)));
%         h{i} = matlabFunction(H{i}(x)/sqrt(ortho_int(@(x)H{i}(x),@(x)H{i}(x),dinfo)));
%     elseif(i==order+2) % Create n+1 order poly for roots
%         temp = matlabFunction(x.*h{i-1}(x));
%         root_poly.H = matlabFunction(x.*h{i-1}(x)-ortho_int(@(x)temp(x),@(x)h{i-1}(x),dinfo)*h{i-1}(x)-...
%             sqrt(ortho_int(@(x)H{i-1}(x),@(x)H{i-1}(x),dinfo))*h{i-2}(x));
%         root_poly.h = matlabFunction(root_poly.H(x)/sqrt(ortho_int(@(x)root_poly.H(x),@(x)root_poly.H(x),dinfo)));
%     else
%         temp = matlabFunction(x.*h{i-1}(x));
%         H{i} = matlabFunction(x.*h{i-1}(x)-ortho_int(@(x)temp(x),@(x)h{i-1}(x),dinfo)*h{i-1}(x)-...
%             sqrt(ortho_int(@(x)H{i-1}(x),@(x)H{i-1}(x),dinfo))*h{i-2}(x));
%         h{i} = matlabFunction(H{i}(x)/sqrt(ortho_int(@(x)H{i}(x),@(x)H{i}(x),dinfo)));
%     end
% end
out.H = H(2:n+1);
out.h = h(2:n+1);
end

function out = ortho_int(p1,p2,dinfo)
syms x   % sig mu
mean = dinfo.mean;
sig = dinfo.std;
dist_type = 2; % 0: normal dist 1: inverse Gaussian dist 2: truncated Gaussian
if(dist_type==0)
    d_func = @(x)1./(sig*sqrt(2*pi))*exp(-((x-mean).^2)/(2*sig^2));
    lim = [-inf inf];
elseif(dist_type==1)
    %alpha = mean^2/sig+2;
    %beta = mean*(mean^2/sig+1);
    %d_func = @(x)(mean.*(mean.^2./sig+1)).^(mean.^2./sig+2)./...
    % gamma(mean.^2./sig+2).*x.^(-(mean.^2./sig+2)-1).*exp(-(mean.*(mean.^2./sig+1)).*x.^(-1));
    d_func = @(x)pdf('InverseGaussian',x,mean,mean^3/sig^2);
    lim = [0 inf];
elseif(dist_type==2)
    if (isfield(dinfo,'trunc')==0)
        error('No truncate info!\n');
    end
    %tmp = @(x)pdf('Normal',x,mean,sig);
    %d_func = truncate(tmp,0,inf);
    d_func = @(x)1./(sig*sqrt(2*pi))*exp(-((x-mean).^2)/(2*sig^2));
    weight = integral(@(x)d_func(x),dinfo.trunc(1),dinfo.trunc(2));
    d_func = @(x)1./(sig*sqrt(2*pi))*exp(-((x-mean).^2)/(2*sig^2))./weight;
    lim = [dinfo.trunc(1) dinfo.trunc(2)];
end
temp = @(x)p1(x).*p2(x).*d_func(x);
%temp = matlabFunction(p1(x)*p2(x)*d_func(x));
%out = integral(@(x)temp(dinfo.mean,dinfo.std,x),lim(1),lim(2));
out = integral(@(x)temp(x),lim(1),lim(2));
end