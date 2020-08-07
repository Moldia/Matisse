%Function: tsneRGB
%2020, Nilsson Lab, Scilifelab
%Objective: Represent in RGB way the expression profile of samples
%Input: Matisse object containing binned data (CELLS/ BINS)
%Output: RGB image our tissue, representing each color a different
%expression pattern.

function tsneRGB(CELLS,varargin)
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




%Filter for 1 variable



% tSNE in MATLAB
% ONLY >=R2018a
seeds = 1e-4*randn(size(cGenes,1), 3);
if ismember(varargin(1),'tsne')
Y = tsne(cGenes, 'NumDimensions', 3, 'NumPCAComponents', 4, 'Perplexity', 30,... %3,4
    'Standardize', 1, 'LearnRate', 1000, 'Verbose', 1, 'InitialY', seeds); 
end
if  ismember(varargin(1),'umap')
   Y=run_umap(cGenes,'n_components',3);
end

if ismember(varargin(1),'pca')
  [coeff, score, latent] = pca(cGenes,'NumComponents',2);
  figure, biplot(coeff(:,1:2), 'Scores', score(:,1:2), 'VarLabels', genes);
  title('top two principle components');
end
 %[uSamples, ~, iSample] = unique(ve);
% [uSamples2, ~, iSample2] = unique(cellfun(@(v) v(1:strfind(v, '_hexbin')-1), table2cell(tableCount(:,1)), 'uni', 0));
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
Yrgb=vertcat(Yrgb,tab)

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
end

VERSE=CELLS.expressionmatrix(vect,:);

CORRELATIONS=zeros(size(Y,2),size(VERSE,2));


for RGB=1:size(Y,2)
for vect=1:size(VERSE,2)
 COS=corrcoef(Y(:,RGB),VERSE(:,vect));
 CORRELATIONS(RGB,vect)=COS(1,2);
end
end    

figure
b=bar(categorical(CELLS.expressionnames),CORRELATIONS')

end