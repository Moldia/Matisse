function [SPOTS]=ISS_getspotsVPRELIMINAR2Y(decoded_file);

data = importdata(decoded_file, ',', 1);
%[name, pos] = getinsitudata(decoded_file);
SPOTS.name=data.textdata(2:end,1);
SPOTS.pos=data.data(:,1:2);



end