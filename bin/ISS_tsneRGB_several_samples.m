% use bin count data to run PCA and tSNE
% Xiaoyan, 2018
function ISS_tsneRGB_several_samples(OUTPUT,SPOTS,minimumexp)


%% do not modify
% import data
tableCount = array2table(OUTPUT.exp);
tableCount.Properties.VariableNames=OUTPUT.genename;
%Filter for 1 column

vect=sum(table2array(tableCount(:,2:end)),2) >= minimumexp;
tableNO=tableCount(sum(table2array(tableCount(:,2:end)),2) < minimumexp, :);
tableCount=tableCount(sum(table2array(tableCount(:,2:end)),2) >= minimumexp, :);





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



%Filter for 1 variable



% tSNE in MATLAB
% ONLY >=R2018a
seeds = 1e-4*randn(size(cGenes,1), 3);
Y = tsne(cGenes, 'NumDimensions', 3, 'NumPCAComponents', 4, 'Perplexity', 40,...
    'Standardize', 1, 'LearnRate', 1000, 'Verbose', 1, 'InitialY', seeds); 
 %[uSamples, ~, iSample] = unique(ve);
% [uSamples2, ~, iSample2] = unique(cellfun(@(v) v(1:strfind(v, '_hexbin')-1), table2cell(tableCount(:,1)), 'uni', 0));
figure, scatter3(Y(:,1),Y(:,2),Y(:,3),10, OUTPUT.sample(vect));
colormap(jet(20));
 title({'tSNE dim reduction to three', 'color-coded by samples'});




% get position
pos = OUTPUT.loc;
pos(:,3)=vect;
pos=array2table(pos);
posno=pos(pos.pos3==0,:);
posno=table2array(posno(:,1:2));
pos=pos(pos.pos3==1,:);
pos=table2array(pos(:,1:2));
% visualize tSNE in RGB (no background)
if ~OUTPUT.hexbinsize;	OUTPUT.hexbinsize = 30;    end
Yrgb = rgbscale(Y);

%Visualization depending on tsne dimension
figure, scatter3(Y(:,1),Y(:,2),Y(:,3),10, Yrgb,'filled');
colormap(jet(20));
title({'tSNE dim reduction to three', 'Color based on axis, according to spatial representation'});



% Visualisation depending on gene expression
% 
for s=1:size(tableCount,2)
genc=table2array(tableCount(:,s));
ti=tableCount.Properties.VariableNames(s);
figure, scatter3(Y(:,1),Y(:,2),Y(:,3),10,genc);
colormap(hsv(100));
title(ti);
end 

%pos=vertcat(pos,posno);
tab=zeros(size(posno,1),3);
tab(:,:)=1;
Yrgb=vertcat(Yrgb,tab)

samplex=unique(OUTPUT.sample);
SUBSAMPLE=OUTPUT.sample(vect);
SUBCOORD=OUTPUT.loc(vect,:);

for element=1:size(samplex)
elename=samplex(element);
DAP=imread(char(OUTPUT.image(element)));
DAP=imresize(DAP,1/cell2mat(OUTPUT.scale(element)));
if size(DAP,3)>1
    DAP=rgb2gray(DAP);
end
figure;
imshow(DAP*5);
hold on
s=scatter(SUBCOORD(ismember(SUBSAMPLE,elename),2),SUBCOORD(ismember(SUBSAMPLE,elename),1),OUTPUT.hexbinsize, Yrgb(ismember(SUBSAMPLE,elename),:),'filled','Marker', 's');
s.MarkerFaceAlpha = 0.70;
hold off
end
%%

% write
csvwrite(fullfile(SPOTS.output_directory, 'tSNE_3D.csv'), Y);
csvwrite(fullfile(SPOTS.output_directory, 'tSNE_initial.csv'), seeds);

end