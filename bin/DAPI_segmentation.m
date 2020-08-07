%Function: DAPI_segmentation
%2020, Nilsson Lab, Scilifelab
%Objective: Define the border of the cells found in a sample based on DAPI
%Input: Matisse object (containing DAPI image)
%Output: Matisse object, that now contains the Cell map and the centroids
%of cells.
%NOTE THAT: SPOTS.segmentation_percent can be modified previously to match
%good segmentation
function [SPOTS]= DAPI_segmentation(SPOTS)

sam=unique(SPOTS.sample);

for eachsam=1:size(sam,1)
sampl=sam(eachsam);

image=SPOTS.images{eachsam};
scale=SPOTS.scale{eachsam};

IDAPI=imread(image);
IDAPI=imresize(IDAPI,1/scale);
if size(IDAPI,3)>1
   IDAPI=rgb2gray(IDAPI); 
end
Dapi = imadjust(IDAPI); % contrast enhancement
ImSz = size(Dapi);
Debug = 1;
%% threshold the map
ThreshVal = prctile(Dapi, SPOTS.segmentation_percent);

bwDapi = imerode(Dapi>ThreshVal, strel('disk', 2));

if Debug
    figure(300)
    subplot(2,1,1);
    imagesc(Dapi); 
    subplot(2,1,2);
    imagesc(bwDapi);
    colormap bone
    fprintf('Threshold = %f\n', ThreshVal);
%    linkaxes();
end
%% find local maxima 
dist = bwdist(~bwDapi);
dist0 = dist;
dist0(dist<5)=0;


ddist = imdilate(dist0, strel('disk', 7));
%clear dist 
impim = imimposemin(-dist0, imregionalmax(ddist));
clear dist0
if Debug
    figure(301);
    subplot(2,1,1)
    imagesc(dist);
    subplot(2,1,2)
    imagesc(impim);
end
%% segment
% remove pixels at watershed boundaries
bwDapi0 = bwDapi;
bwDapi0(watershed(impim)==0)=0;

% assign all pixels a label
labels = uint32(bwlabel(bwDapi0));
[d, idx] = bwdist(bwDapi0);

% now expand the regions by a margin
CellMap0 = zeros(ImSz, 'uint32');
Expansions = (d<10);
CellMap0(Expansions) = labels(idx(Expansions));

% get rid of cells that are too small
rProps0 = regionprops(CellMap0); % returns XY coordinate and area
BigEnough = [rProps0.Area]>200;
NewNumber = zeros(length(rProps0),1);
NewNumber(~BigEnough) = 0;
NewNumber(BigEnough) = 1:sum(BigEnough);
CellMap = CellMap0;
CellMap(CellMap0>0) = NewNumber(CellMap0(CellMap0>0));


U2=CellMap;
 blobs = unique(U2);
 blobs = blobs(blobs~=0);
 size(blobs);
XY=fliplr(vertcat(regionprops(CellMap).Centroid));

SPOTS.centroidpos{eachsam}=XY;
SPOTS.centroidID{eachsam}=blobs;
SPOTS.Cellmaps{eachsam}=CellMap;
end

end

