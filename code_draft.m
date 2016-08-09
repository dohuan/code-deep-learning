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