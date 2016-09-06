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






