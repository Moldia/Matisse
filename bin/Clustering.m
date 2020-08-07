%Function: Clustering
%2020, Nilsson Lab, Scilifelab
%Objective: Cluster bins/cells data for finding cell type or expression
%patterns
%Input: Matisse object + Number of clusters + method 
  %Method can have two options: 'hierarchical' for hierarchical clustering
  %or 'kmeans' for kmeans clustering
%Output: Matisse object given as an input , but clustered . 
function[CELLS] = Clustering(CELLS,n_clust,METHOD)

%% do not modify
% import data
tableCount = array2table(CELLS.expressionmatrix);
tableCount.Properties.VariableNames=CELLS.expressionnames;

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

if strcmp(METHOD,'kmeans')
% k-means clustering
disp('Starting kmeans clustering with 500 replicates..');
[iCluster, centroid] = kmeans(cGenes, n_clust,...
    'Distance', 'sqeuclidean', 'Replicates', 500);
end

if strcmp(METHOD,'hierarchical')
Z = linkage(cGenes,'ward');
iCluster = cluster(Z,'Maxclust',n_clust)';
end

if strcmp(METHOD,'umap')
lastwarn('Keep in mind that umap clustering does not consider number of clusters')
[CELLS.UMAP_clustering,COORD,iCluster]=run_umap(cGenes);
figure
gscatter(CELLS.UMAP_clustering(:,1),CELLS.UMAP_clustering(:,2),iCluster);
title('UMAP representation of clusters');
xlabel('X UMAP 1')
ylabel('Y UMAP 2')
figure
gscatter(CELLS.UMAP_clustering(:,1),CELLS.UMAP_clustering(:,2),CELLS.sample);
title('UMAP by samples');
xlabel('X UMAP 1')
ylabel('Y UMAP 2')
end
   
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

% figure
% bh = barh(centroid');
% set(gca, 'ytick', 1:nnz(isSelected), 'yticklabel', cNames(isSelected),...
%     'ylim',[0 nnz(isSelected)+1], 'fontsize', 5, 'ydir', 'reverse');
% box off
% xlabel('bin count')
% title('centroid location of clusters')
% legend(catstrnum('Cluster ', 1:max(iCluster)))

%Now we display images itself

samples=unique(CELLS.sample);
ALLCELLS=struct();

for sam=1:size(samples,1)
numb=samples(sam);
goodsamples=ismember(CELLS.sample,numb);
if size(samples,1)>1;
im=imread(char(CELLS.images(sam)));
else
im=imread(char(CELLS.images(sam)));
end
ima=imresize(im,1/cell2mat(CELLS.scale(sam)));
figure
if size(ima,3)>1
imshow(rgb2gray(ima*10));
else 
imshow(ima*10);  
end    
hold on
g=gscatter(CELLS.location(goodsamples,1),CELLS.location(goodsamples,2),iCluster(goodsamples),[],[],10);
hold off
end

CELLS.Clustering_ID=iCluster;
CELLS.spotname=cellstr(num2str(iCluster'));

end 
