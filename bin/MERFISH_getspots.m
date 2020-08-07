function [SPOTS]=MERFISH_getspots(decoded_file,barcodes);

data = importdata(decoded_file, ',', 1);
barcod = importdata(barcodes, ',', 1);
%[name, pos] = getinsitudata(decoded_file);
data.data(:,2)=data.data(:,2)-min(data.data(:,2));
data.data(:,3)=data.data(:,3)-min(data.data(:,3));
SPOTS.pos=data.data(:,2:3);
barco={};
for elem=1:size(data.data,1)
 barco(elem)=barcod.textdata(find(barcod.data(:,1)==data.data(elem,1))+1,1);  
end
SPOTS.name=barco';
end