%Function: show_spots
%2020, Nilsson Lab, Scilifelab
%Objective: Plot the gene reads/spots/cells for visualization purposes
%Input: Matisse object.
%Output: One matlab figure for each sample found in the Matisse object.

function  show_spots(SPOTS,varargin)
samples=unique(SPOTS.sample);
ALLCELLS=struct();

for sam=1:size(samples,1)
numb=samples(sam);
goodsamples=ismember(SPOTS.sample,numb);
figure
if nargin>1
if ~ismember(varargin{1},'blank')
imshow(imresize(imread(SPOTS.images{sam}),1/SPOTS.scale{sam})*10);
hold on
end
else 
imshow(imresize(imread(SPOTS.images{sam}),1/SPOTS.scale{sam})*10);
hold on
end

gscatter(SPOTS.location(goodsamples,1),SPOTS.location(goodsamples,2),SPOTS.spotname(goodsamples,:),[],'o*xdsph');
hold off
end
end