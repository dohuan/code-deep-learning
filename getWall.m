function wall=getWall(ID,scan)
filePath = 'C:\Users\dohuan\Documents\GitHub\wall-thickness\data\';
option = struct('begslice',2,'endslice',4,'L',8,'x',7,'N',20,'pts',36,'spacing',1);
directory.img = [filePath ID ' ' scan ' AAA'];
directory.mask = [filePath ID ' ' scan ' AAA Rough Lumen Mask'];

data = featureExtract(directory,option);


end