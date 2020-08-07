%Function: LowDimensionalRGB
%2020, Nilsson Lab, Scilifelab
%Objective: Represent in RGB way the expression profile of samples
%Input: Matisse object containing binned data (CELLS/ BINS)
%Output: Low dimensional image our tissue, representing each color a different
%expression pattern.

function LowDimensionalRGB(CELLS,varargin)
%% do not modify
% import data
tableCount = array2table(CELLS.expressionmatrix);
tableCount.Properties.VariableNames=CELLS.expressionnames;
%Filter
%r for 1 column
vect=sum(table2array(tableCount(:,2:end)),2) >= CELLS.count_threshold;
tableNO=tableCount(sum(table2array(tableCount(:,2:end)),2) < CELLS.count_threshold, :);
tableCount=tableCount(sum(table2array(tableCount(:,2:end)),2) >= CELLS.count_threshold, :);





% original column names
cNames = tableCount.Properties.VariableNames;

% make sure there are no multiple entries of the same gene
assert(numel(cNames)== numel(unique(cNames)),...
    'Column names are not unique!')

% create checkboxes and get selected values
cbValues = checkboxes(cNames);
idx = cellfun(@(v) find(strcmp(v, cNames)), cbValues(:,1));
isSelected = false(numel(idx), 1);
isSelected(idx) = cell2mat(cbValues(:,2));
cGenes = table2array(tableCount(:,isSelected));
genes = cNames(isSelected)'

% PCA
[coeff, score, latent] = pca(cGenes);
% visualize first two components
figure, biplot(coeff(:,1:2), 'Scores', score(:,1:2), 'VarLabels', genes);
title('top two principle components');


% Three dimensional reductions applied in MATLAB
% ONLY >=R2018a
seeds = 1e-4*randn(size(cGenes,1), 3);
if ismember(varargin(1),'tsne')
Y = tsne(cGenes, 'NumDimensions', 3, 'NumPCAComponents', 4, 'Perplexity', 30,... %3,4
    'Standardize', 1, 'LearnRate', 1000, 'Verbose', 1, 'InitialY', seeds); 
elseif  ismember(varargin(1),'umap')
   Y=run_umap(cGenes,'n_components',3);
elseif ismember(varargin(1),'pca')
  [coeff, score, latent] = pca(cGenes,'NumComponents',3);
  figure, biplot(coeff(:,1:2), 'Scores', score(:,1:2), 'VarLabels', genes);
  title('top two principle components');
  Y=score(:,1:3);
else
    error('Error:Please specify a correct low dimensional method');
end

figure, scatter3(Y(:,1),Y(:,2),Y(:,3),10, categorical(CELLS.sample(vect)));
colormap(jet(20));
 title({'tSNE dim reduction to three', 'color-coded by samples'});




% get position
pos = CELLS.location;
pos(:,3)=vect;
pos=array2table(pos);
posno=pos(pos.pos3==0,:);
posno=table2array(posno(:,1:2));
pos=pos(pos.pos3==1,:);
pos=table2array(pos(:,1:2));
% visualize tSNE in RGB (no background)
Yrgb = rgbscale(Y);

%Visualization depending on tsne dimension
figure, scatter3(Y(:,1),Y(:,2),Y(:,3),10, Yrgb,'filled');
colormap(jet(20));
title({'tSNE dim reduction to three', 'Color based on axis, according to spatial representation'});


% Visualisation depending on gene expression
if ismember(varargin(2),'genexpression')
for s=1:size(tableCount,2)
genc=table2array(tableCount(:,s));
ti=tableCount.Properties.VariableNames(s);
figure, scatter3(Y(:,1),Y(:,2),Y(:,3),10,genc);
colormap(hsv(100));
title(ti);
end 
end

%pos=vertcat(pos,posno);
tab=zeros(size(posno,1),3);
tab(:,:)=1;
Yrgb=vertcat(Yrgb,tab);

samplex=unique(CELLS.sample);
SUBSAMPLE=CELLS.sample(vect);
SUBCOORD=CELLS.location(vect,:);

for element=1:size(samplex,1)
elename=samplex(element);
DAP=imread(CELLS.images{element});
DAP=imresize(DAP,1/cell2mat(CELLS.scale(element)));

if size(DAP,3)>1
   DAP=rgb2gray(DAP);
end
figure;
imshow(DAP*5);
hold on
s=scatter(SUBCOORD(ismember(SUBSAMPLE,elename),1),SUBCOORD(ismember(SUBSAMPLE,elename),2),CELLS.hexbin_size, Yrgb(ismember(SUBSAMPLE,elename),:),'filled','Marker', 's');
s.MarkerFaceAlpha = 0.70;
hold off
%%MODIFICATION OF THE IMAGE IN ORDER TO VISUALIZE IT NICELY
rgbImage = cat(3, zeros(size(DAP)),zeros(size(DAP)),zeros(size(DAP)));
%ismember(SUBSAMPLE,elename)
PLOTCOORD=SUBCOORD(ismember(SUBSAMPLE,elename),:);
PLOTVALUES=Yrgb(ismember(SUBSAMPLE,elename),:);
PLOTVALUES=PLOTVALUES;
INTERV=max(round((PLOTCOORD(2,1)-PLOTCOORD(1,1))),round((PLOTCOORD(2,2)-PLOTCOORD(1,2))));
%%%%%%%%%%%%%%THIS PART PLOTS THE REAL IMAGE%%%%%%%%%%%%%%%%%%%%
for elis=1:size(PLOTCOORD,1)
disp([num2str(elis),'/',num2str(size(PLOTCOORD,1))]);
for x=round(PLOTCOORD(elis,1)-INTERV):round(PLOTCOORD(elis,1)+INTERV)
 for y=round(PLOTCOORD(elis,2)-INTERV):round(PLOTCOORD(elis,2)+INTERV)
if x>0 & x<size(DAP,2) & y>0 & y<size(DAP,1);
% disp(['---------',num2str(x)]);
rgbImage(y,x,:)=PLOTVALUES(elis,:);
end
end
end
end
figure;
l=imshow(rgbImage);
%  h.alpha = DAP>0;
 hold on
h=imshow(DAP*5);
set(h, 'AlphaData', 0.2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


VERSE=CELLS.expressionmatrix(vect,:);
CORRELATIONS=zeros(size(Y,2),size(VERSE,2));

%We create the bar plot showing the weight of each color in the 3
%dimensional representation

for RGB=1:size(Y,2)
for vect=1:size(VERSE,2)
 COS=corrcoef(Y(:,RGB),VERSE(:,vect));
 CORRELATIONS(RGB,vect)=COS(1,2);
end
end    

figure
b=bar(categorical(CELLS.expressionnames),CORRELATIONS','FaceColor','flat');
b(1).FaceColor = [1 0 0]
b(2).FaceColor = [0 0.5 0]
b(3).FaceColor = [0 0 0.7]
end