function [SPOTS]=osmFISH_getspots(decoded_file);

data = importdata(decoded_file, ',', 1);
SPOTS.pos=data.data(:,3:4);
SPOTS.name=data.textdata(2:end,3);

end