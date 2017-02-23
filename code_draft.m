% --- prepare data from DCM files
dcmfiles = dir('./data/P11 S04 AAA/*.dcm');
data = [];
for i=1:length(dcmfiles)
	fprintf([dcmfiles(i).name '\n']);
	img = int16(dicomread(['./data/P11 S04 AAA/' dcmfiles(i).name]));
	info = dicominfo(['./data/P11 S04 AAA/' dcmfiles(i).name]);
	img = double(img*info.RescaleSlope + info.RescaleIntercept);
	img = (img-min(img(:)))./(max(img(:))-min(img(:)));
	img = imcrop(img,[121 109 300 300]);
	img = img(1:3:end,1:3:end);
	imshow(img)
	data = [data; reshape(img,1,[])];
end
save('data','data');


% --- Check stats of scan times
time = [];
for i=1:14
	if (length(Time{i})>2)
		time = [time; Time{i}(end-2:end)'];
	end
end
time = time.*30;


% --- Generate UNIFORMLY distributed parameters k, sigma, and mu
% run ONLY once then distribute value to different PCs
range.k = [0.196 0.287];
range.sigma = [1.72 2.70];
range.mu = [8.77 9.71];
ran_k = range.k(1)+(range.k(2)-range.k(1)).*rand(8,1);
ran_sigma = range.sigma(1)+(range.sigma(2)-range.sigma(1)).*rand(8,1);
ran_mu = range.mu(1)+(range.mu(2)-range.mu(1)).*rand(8,1);
[k sig mu] = meshgrid(ran_k,ran_sigma,ran_mu);
param = [k(:) sig(:) mu(:)];



% --- Merge GR data together
% 1. load data from file
% 2. subsample by a time period of 12 months
% NOTE: save data as SINGLE (4 bytes) to save space
fileList = dir('./data/GR_H_081716/*.mat');
data.train_x = [];
data.train_y  = [];
fnum = 3; % take [fnum] scans in the past as features, predict the scan [fnum+1]-th
for i=1:length(fileList)
	loadFile = ['./data/GR_H_081716/' fileList(i).name];
	fprintf(['Loading data from: ' loadFile '\n']);
	load(loadFile);
	t_offset = 12*(fnum+1)*30/5; % t_offset = 12(months)*(fnum+1)(scans)*30(days)/5(days)
	ran_ix = randi([1,size(result,1)-t_offset],100,1); 
	t_step = 12*30/5; % =72
	for j=1:length(ran_ix)
		tmp_x = [];
		stoken = (rand()>.7); % shift token
		for k=1:fnum+1
			ix = ran_ix(j)+t_step*(k-1);
			if (stoken ==1 )
				result_tmp = (curve_shifter(single(result(ix,:)'), .7))';
			else
				result_tmp = single(result(ix,:));
			end
			if k~=fnum+1
				tmp_x = [tmp_x, result_tmp];
			else
				
				tmp_y = result_tmp;
			end
		end
		data.train_x = [data.train_x;tmp_x];
		data.train_y = [data.train_y; tmp_y];
	end
end
data.fnum = fnum;




% --- Load and NORMALIZED REAL data 
fileList = dir('./data/real_scans/*.mat');
dataout.ft_x   = []; % ft = 'fine tune'
dataout.ft_y   = [];
dataout.test_x = [];
dataout.test_y = [];
fnum = 3; % take [fnum] scans in the past as features, predict the scan [fnum+1]-th
% NOTE: make sure fnum is the SAME as GR data
dataout.fnum = fnum;
count = 1;
for i = 1:length(fileList)
	loadFile = ['./data/real_scans/' fileList(i).name];
	load(loadFile);
	if (length(data)<fnum+1)
		fprintf('Discarded: insufficient data\n');
	else
		grp_num = (length(data)-(fnum+1)+1);
		fprintf('Possible data size: %d\n',grp_num);
		grp_ix = 1:length(data);
		grp_ix = grp_ix';
		for j=1:grp_num
			ix_tmp = circshift(grp_ix,-j+1);
			data_grp = data(ix_tmp(1:fnum+1));
			tmp_x = [];
			for k=1:length(data_grp)
				p = polyfit(data_grp(k).pos, data_grp(k).maxd,7);
				x = linspace(min(data_grp(k).pos), max(data_grp(k).pos), 221);
				tmp = polyval(p,x);
				if (k~=length(data_grp))
					tmp_x = [tmp_x, tmp];
				else
					tmp_y = tmp;
				end
				
			end
			% --- Un-comment following two lines if normalized
			%tmp_x = tmp_x./max(tmp_x);
			%tmp_y = tmp_y./max(tmp_y);
			if(j~=grp_num)
				dataout.ft_x = [dataout.ft_x; tmp_x];
				dataout.ft_y = [dataout.ft_y; tmp_y];
			else
				dataout.test_x = [dataout.test_x; tmp_x];
				dataout.test_y = [dataout.test_y; tmp_y];
			end
		end
	end
end
clear data
data = dataout;


% --- Load and NORMALIZED REAL data for CROSS VALIDATION
fileList = dir('./data/real_scans/*.mat');
fnum = 3; % take [fnum] scans in the past as features, predict the scan [fnum+1]-th
% NOTE: make sure fnum is the SAME as GR data
count = 1;
for i = 1:length(fileList)
	loadFile = ['./data/real_scans/' fileList(i).name];
	load(loadFile);
	if (length(data)<fnum+1)
		fprintf('Discarded: insufficient data\n');
	else
		
		dataout(count).ft_x   = []; % ft = 'fine tune'
		dataout(count).ft_y   = [];
		dataout(count).test_x = [];
		dataout(count).test_y = [];
		dataout(count).fnum = fnum;
		
		grp_num = (length(data)-(fnum+1)+1);
		fprintf('Possible data size: %d\n',grp_num);
		grp_ix = 1:length(data);
		grp_ix = grp_ix';
		for j=1:grp_num
			ix_tmp = circshift(grp_ix,-j+1);
			data_grp = data(ix_tmp(1:fnum+1));
			tmp_x = [];
			for k=1:length(data_grp)
				p = polyfit(data_grp(k).pos, data_grp(k).maxd,7);
				x = linspace(min(data_grp(k).pos), max(data_grp(k).pos), 221);
				tmp = polyval(p,x);
				if (k~=length(data_grp))
					tmp_x = [tmp_x, tmp];
				else
					tmp_y = tmp;
				end
				
			end
			% --- Un-comment following two lines if normalized
			%tmp_x = tmp_x./max(tmp_x);
			%tmp_y = tmp_y./max(tmp_y);
			if(j~=grp_num)
				dataout(count).ft_x = [dataout(count).ft_x; tmp_x];
				dataout(count).ft_y = [dataout(count).ft_y; tmp_y];
			else
				dataout(count).test_x = [dataout(count).test_x; tmp_x];
				dataout(count).test_y = [dataout(count).test_y; tmp_y];
			end
		end
		count = count + 1;
	end
end
clear data
data = dataout;






A = 1:10;
A = A';
B(:,1) = A;
for i=2:length(A)
	%B(:,i) = circshift(B(:,i-1),-1);
	B(:,i) = circshift(A,-i+1);
end







% --- Test polyfit with patient H data
% Assume already load patient H data
index = 4;
p = polyfit(H(index).pos, H(index).maxd,7);
x = linspace(min(H(index).pos), max(H(index).pos), 221);
y = polyval(p,x);
figure(1)
hold on
plot(H(index).pos, H(index).maxd);
plot(x,y);
hold off


% --- Fix patient K data NaN error
for i=1:length(data)
	h = isnan(data(i).maxd);
	ix = find(h==1);
	data(i).maxd(ix) = [];
	data(i).pos(ix) = [];
end


% --- Find out the universal maximum diameter
load ./data/dataGR-augmented
uni_max = 0;
for i=1:size(data.train_y,1)
	if (uni_max<max(data.train_y(i,:)))
		uni_max = max(data.train_y(i,:));
	end
end
fprintf('Max dia of GR data: %.3f\n',uni_max);
load ./data/dataREAL_unnormalized
uni_max = 0;
for i=1:size(data.test_y,1)
	if (uni_max<max(data.test_y(i,:)))
		uni_max = max(data.test_y(i,:));
	end
end
fprintf('Max dia of REAL data: %.3f\n',uni_max);



% --- Scale all data by the universal max dia (found by the code above)
uni_max = 7.1;
load ./data/dataGR-augmented
train_x = data.train_x./uni_max;
train_y = data.train_y./uni_max;
save('./data/data_train','train_x','train_y','uni_max');

load ./data/dataREAL_unnormalized
data.ft_x = data.ft_x./uni_max;
data.ft_y = data.ft_y./uni_max;
data.test_x = data.test_x./uni_max;
data.test_y = data.test_y./uni_max;
save('./data/dataREAL_normalized_uni_max','data')




h = repmat(rbm.c', opts.batchsize, 1) + v1 * rbm.W';


% --- plot histogram of GR and CM side by side
GR = load('./data/dataGR-augmented');
CM = load('./data/dataCM-augmented');
gr = GR.data.train_y;
cm = CM.data.train_y;
hold on
for i=1:size(gr,2)
	[h,x] = hist(gr(:,i));
	[h1,x1] = hist(cm(:,i));
	h=h./(max(h));
	h1=h1./max(h1);
	plot3(x,i*ones(1,length(x)),h,'b-');
	plot3(x1,i*ones(1,length(x1)),h1,'r--');
end
hold off
xlabel('Maximal diameter')
ylabel('Spatial site')
zlabel('Probability')



% --- plot epoch
epoch = [1 10 100 150 200 250 300 500 700 900];
RMSE = [7.63 6.24 4.83 5.07 4.81 5.03 5.07 4.68 4.69 4.75];
ttime = [29.99 29.58 34.58 36.38 38.24 41.20 42.92 51.38 64.23 69.44];
[Ax,h1,h2] = plotyy(epoch,RMSE,epoch,ttime,'plot');
h1.LineStyle = '-';
h2.LineStyle = '--';
h1.LineWidth = 2;
h2.LineWidth = 2;
ylabel(Ax(1),'RMSE');
ylabel(Ax(2),'Training time (seconds)');
xlabel('Number of epochs');
set(gca,'FontSize',16);





% --- Compare proposed method with linear mixed-effects
ME = load('mvgress_try');
%DL = load('results020117');
DL = load('results_dropout');
true = DL.test_y;
estME = ME.est_y;
estDL = DL.est;
scaleDL = DL.scale_test;

patientID{1} = 'P1'; %G
patientID{2} = 'P2'; %H
patientID{3} = 'P3'; %J
patientID{4} = 'P4'; %K
patientID{5} = 'P5'; %P11
patientID{6} = 'P6'; %P12
patientID{7} = 'P7'; %P13


for i=1:7
	fprintf('RMSE of patient %s: %.2f (DL) %.2f (ME)\n',patientID{i},rmseCal(estDL(i,:).*scaleDL(i),true(i,:).*DL.scaleTrackTest(i,end)),rmseCal(estME(i,:),true(i,:).*DL.scaleTrackTest(i,end)));
end

for i=1:7
	fprintf('RMSE of patient %s (normalized): %.2f (DL) %.2f (ME)\n',patientID{i},rmseCal(estDL(i,:),true(i,:)),rmseCal(estME(i,:)./max(estME(i,:)),true(i,:)));
end

count = 1;
figure(2)
for i=[1 2 3 5 6 7]
    subplot(2,3,count)
    hold on
    estPlot = estDL(i,:).*scaleDL(i);
    %estPlot = smooth(estPlot,.1,'lowess');
    %plot(est(i,:).*scale_test(i),'g-.','LineWidth',2)
    plot(true(i,:)*DL.scaleTrackTest(i,end),'k-','LineWidth',2)
    plot(estPlot,'b--','LineWidth',2);
	plot(estME(i,:),'r-.','LineWidth',2)
    hold off
    legend('true','DL prediction','ME prediction')
    title(['P' num2str(count)])
    box on
    axis tight
	count = count + 1;
end
%set(gca,'FontSize',16);



% --- Compare proposed method (DROPOUT and NONDROPOUT) with linear mixed-effects
ME = load('mvgress_try');
DL = load('results020117_published');
DO = load('results_dropout');
true = DL.test_y;
estME = ME.est_y;
estDL = DL.est;
scaleDL = DL.scale_test;
estDO = DO.est;
scaleDO = DO.scale_test;


patientID{1} = 'P1'; %G
patientID{2} = 'P2'; %H
patientID{3} = 'P3'; %J
patientID{4} = 'P4'; %K
patientID{5} = 'P5'; %P11
patientID{6} = 'P6'; %P12
patientID{7} = 'P7'; %P13


for i=1:7
	fprintf('RMSE of patient %s: %.2f (DL) %.2f (DO) %.2f (ME)\n',patientID{i},rmseCal(estDL(i,:).*scaleDL(i),true(i,:).*DL.scaleTrackTest(i,end)),rmseCal(estDO(i,:).*scaleDO(i),true(i,:).*DL.scaleTrackTest(i,end)),rmseCal(estME(i,:),true(i,:).*DL.scaleTrackTest(i,end)));
end

for i=1:7
	fprintf('RMSE of patient %s (normalized): %.2f (DL) %.2f (DO) %.2f (ME)\n',patientID{i},rmseCal(estDL(i,:),true(i,:)),rmseCal(estDO(i,:),true(i,:)),rmseCal(estME(i,:)./max(estME(i,:)),true(i,:)));
end

count = 1;
figure(2)
for i=[1 2 3 5 6 7]
    subplot(2,3,count)
    hold on
    estPlot = estDL(i,:).*scaleDL(i);
	estPlotDO = estDO(i,:).*scaleDO(i);
    %estPlot = smooth(estPlot,.1,'lowess');
    %plot(est(i,:).*scale_test(i),'g-.','LineWidth',2)
    plot(true(i,:)*DL.scaleTrackTest(i,end),'k-','LineWidth',2)
    plot(estPlot,'b--','LineWidth',2);
	plot(estPlotDO,'c:','LineWidth',2);
	plot(estME(i,:),'r-.','LineWidth',2)
    hold off
    legend('true','DL prediction','DO prediction','ME prediction')
    title(['P' num2str(count)])
    box on
    axis tight
	count = count + 1;
end
%set(gca,'FontSize',16);





% ----
h= rawdata(ix1,:);
h1 = (curve_shifter(h',.7))';
h2 = curve_merger(h,rawdata(ix2,:),1,ratio);
hold on
plot(h,'b','LineWidth',2);
plot(h1,'r--','LineWidth',2);
plot(h2,'k:','LineWidth',2);
box on
axis tight
xlabel('centerline (spatial site unit)');
ylabel('maximal diameter (cm)');
set(gca,'FontSize',16);




% --- plot the weights
load results020117_published

figure(1)
subplot(2,1,1)
title('RBM (pre-trained)')
hold on
for i=1:size(dbn.rbm{1}.W,1)
plot(dbn.rbm{1}.W(i,:));
end
hold off
box on
axis tight
ylabel('weights')
xlabel('layer 1 units')
subplot(2,1,2)
title('NN (fine-tuned)')
hold on
for i=1:size(nn.W{1},1)
plot(nn.W{1}(i,:));
end
ylabel('weights')
xlabel('layer 1 units')
hold off
box on
axis tight

figure(2)
subplot(2,1,1)
title('RBM (pre-trained)')
hold on
for i=1:size(dbn.rbm{2}.W,1)
plot(dbn.rbm{2}.W(i,:));
end
hold off
box on
axis tight
ylabel('weights')
xlabel('layer 2 units')
subplot(2,1,2)
title('NN (fine-tuned)')
hold on
for i=1:size(nn.W{2},1)
plot(nn.W{2}(i,:));
end
ylabel('weights')
xlabel('layer 2 units')
hold off
box on
axis tight

figure(3)
subplot(2,1,1)
title('RBM (pre-trained)')
hold on
for i=1:size(dbn.rbm{3}.W,1)
plot(dbn.rbm{3}.W(i,:));
end
hold off
box on
axis tight
ylabel('weights')
xlabel('layer 3 units')
subplot(2,1,2)
title('NN (fine-tuned)')
hold on
for i=1:size(nn.W{3},1)
plot(nn.W{3}(i,:));
end
ylabel('weights')
xlabel('layer 3 units')
hold off
box on
axis tight




% -----------

subplot(1,4,[1 2 3])
title('Training data')
plot(data.train_x(1,:),'LineWidth',2)
axis tight
ylabel('Maximum diameter (cm)')
ylim([2.0 2.8])
set(gca,'FontSize',16)
subplot(1,4,4)
title('Testing data')
plot(data.train_y(1,:),'LineWidth',2)
axis tight
ylabel('Maximum diameter (cm)')
set(gca,'FontSize',16)


% --- plot results for CV
load ./results_CV
patientID{1} = 'P1'; %G
patientID{2} = 'P2'; %H
patientID{3} = 'P3'; %J
patientID{4} = 'P4'; %K
patientID{5} = 'P5'; %P11
patientID{6} = 'P6'; %P12
patientID{7} = 'P7'; %P13

count = 1;
figure(1)
for i=[1 2 3 5 6 7]
    subplot(2,3,count)
    hold on
    plot(mdc_est(i,:),'b--','LineWidth',2)
    plot(mdc_true(i,:),'r','LineWidth',2);
    hold off
    legend('Predicted','true')
    title(['P' num2str(count)])
    box on
    axis tight
	count = count + 1;
end