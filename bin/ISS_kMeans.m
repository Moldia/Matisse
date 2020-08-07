% use count data to run k-means clustering
% select which columns to include
% Xiaoyan, 2017
function[ALLCELLS] = ISS_kMeans(ALLEXPRESSION,num_clusters)

%% do not modify
% import data
tableCount = array2table(ALLEXPRESSION.exp);
tableCount.Properties.VariableNames=ALLEXPRESSION.genename;

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

% k-means clustering
disp('Starting kmeans clustering with 500 replicates..');
[iCluster, centroid] = kmeans(cGenes, num_clusters,...
    'Distance', 'sqeuclidean', 'Replicates', 500);

% heatmap of normalized counts and barplot of cluster centroid position
figure;

[~, idxSort] = sort(iCluster);
[~, idxFirst] = unique(iCluster(idxSort));
imagesc(cGenes(idxSort,:)');
hold on;
plot(repmat(idxFirst', 2, 1), repmat([0; nnz(isSelected)+1], 1, numel(idxFirst)),...
    'black', 'linewidth', 1);
xlabel('bin');
set(gca, 'ytick', 1:nnz(isSelected), 'yticklabel', cNames(isSelected), 'fontsize', 8);
title('bin count data')
colorbar

figure
bh = barh(centroid');
set(gca, 'ytick', 1:nnz(isSelected), 'yticklabel', cNames(isSelected),...
    'ylim',[0 nnz(isSelected)+1], 'fontsize', 5, 'ydir', 'reverse');
box off
xlabel('bin count')
title('centroid location of clusters')
legend(catstrnum('Cluster ', 1:max(iCluster)))

%Now we display images itself

samples=unique(ALLEXPRESSION.sample);
ALLCELLS=struct();

for sam=1:size(samples,1)
numb=samples(sam);
goodsamples=ismember(ALLEXPRESSION.sample,numb);
im=imread(char(ALLEXPRESSION.image(sam)));
ima=imresize(im,1/cell2mat(ALLEXPRESSION.scale(sam)));
figure
imshow(rgb2gray(ima*10));
hold on
g=gscatter(ALLEXPRESSION.loc(goodsamples,2),ALLEXPRESSION.loc(goodsamples,1),iCluster(goodsamples),[],[],10);

CELLS.matrixcount=tableCount(goodsamples,:);
CELLS.geneames=ALLEXPRESSION.genename;
CELLS.name=iCluster(goodsamples)';
CELLS.pos=[ALLEXPRESSION.loc(goodsamples,2),ALLEXPRESSION.loc(goodsamples,1)];
CELLS.image=ALLEXPRESSION.image(sam);
CELLS.scale=ALLEXPRESSION.scale(sam);
%ALLCELLS{sam}=CELLS;
ALLCELLS = setfield(ALLCELLS,['sample',num2str(sam)],CELLS);
end



end 
