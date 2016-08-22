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
fnum = 4; % take [fnum] scans in the past as features, predict the scan [fnum+1]-th
for i=1:length(fileList)
	loadFile = ['./data/GR_H_081716/' fileList(i).name];
	fprintf(['Loading data from: ' loadFile '\n']);
	load(loadFile);
	t_offset = 12*(fnum+1)*30/5; % t_offset = 12(months)*(fnum+1)(scans)*30(days)/5(days)
	ran_ix = randi([1,size(result,1)-t_offset],100,1); 
	t_step = 12(months)*30(days)/5(days); % =72
	for j=1:length(ran_ix)
		tmp_x = [];
		for k=1:fnum+1
			ix = ran_ix(j)+t_step*(k-1);
			if k~=fnum+1
				tmp_x = [tmp_x, single(result(ix,:))];
			else
				tmp_y = single(result(ix,:));
			end
		end
		data.train_x = [data.train_x;tmp_x];
		data.train_y = [data.train_y; tmp_y];
	end
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
















