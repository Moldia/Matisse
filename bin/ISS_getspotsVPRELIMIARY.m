function [SPOTS]=ISS_getspotsVPRELIMINARY2(decoded_file);

data = importdata(decoded_file, ',', 1);
%[name, pos] = getinsitudata(decoded_file);
SPOTS.name=data.textdata(2:end,2);
SPOTS.pos=data.data(:,1:2);



end