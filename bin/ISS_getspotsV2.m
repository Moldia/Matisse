function [SPOTS]=ISS_getspots(decoded_file,varargin);

if nargin <1
%data = importdata(decoded_file, ',', 1);
data=readtable(decoded_file);
%[name, pos] = getinsitudata(decoded_file);
SPOTS.name=data(2:end,1);
SPOTS.pos=data(:,1:2);
else
%data = importdata(decoded_file, ',', 1);
data=readtable(decoded_file);
%[name, pos] = getinsitudata(decoded_file);
SPOTS.name=data(:,varargin{2});
SPOTS.name=table2cell(SPOTS.name);
SPOTS.pos=data(:,varargin{4});   
SPOTS.pos=table2array(SPOTS.pos);
end

end