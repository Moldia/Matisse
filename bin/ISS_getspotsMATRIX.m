function [SPOTS]=ISS_getspotsMATRIX(decoded_file);

data = importdata(decoded_file, ',', 1);
%[name, pos] = getinsitudata(decoded_file);
SPOTS.name=data.textdata(2:end,1);
SPOTS.pos=data.data;



end