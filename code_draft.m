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
data.train = [];
data.test  = [];
for i=1:length(fileList)
	loadFile = ['./data/GR_H_081716/' fileList(i).name];
	fprintf(['Loading data from: ' loadFile '\n']);
	load(loadFile);
	ran_ix = randi([1,size(result,1)-288],100,1); % 288 = 12(months)*4(scans)*30(days)/5(days)
	for j=1:length(ran_ix)
		tmp_train = [];
		for k=1:5
			ix = ran_ix(j)+72*(k-1); % 72 = 12(months)*30(days)/5(days)
			if k~=5
				tmp_train = [tmp_train, single(result(ix,:))];
			else
				tmp_test = single(result(ix,:));
			end
		end
		data.train = [data.train;tmp_train];
		data.test = [data.test; tmp_test];
	end
end


















