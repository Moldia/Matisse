%Function: MergeObjects
%2020, Nilsson Lab, Scilifelab
%Objective: Combine two matisse objects in a single one
%Input: Two different Matisse objects
%Output: Combined Reference Mastisse object
function [FIN] = MergeObjects(SPOTS,CORTEX_SPOTS)
FIN=SPOTS;
FIN.location=[FIN.location;CORTEX_SPOTS.location];
FIN.spotname=[FIN.spotname;CORTEX_SPOTS.spotname];
FIN.sample=[FIN.sample;CORTEX_SPOTS.sample];
FIN.images=[FIN.images,CORTEX_SPOTS.images];
FIN.scale=[FIN.scale,CORTEX_SPOTS.scale];
FIN.expressionmatrix=vertcat(FIN.expressionmatrix,CORTEX_SPOTS.expressionmatrix);

end
