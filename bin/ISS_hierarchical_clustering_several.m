function [ALLCELLS]=ISS_hierarchical_clustering_several(ALLEXPRESSION,SPOTS,n_clust)
matCounts=ALLEXPRESSION.exp;
Z = linkage(matCounts(:,:),'ward');
c = cluster(Z,'Maxclust',n_clust);
%figure
%X=tsne(matCounts, 'NumDimensions', 3, 'NumPCAComponents', 4);
%scatter3(X(:,1),X(:,2),X(:,3),10,c)
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
g=gscatter(ALLEXPRESSION.loc(goodsamples,2),ALLEXPRESSION.loc(goodsamples,1),c(goodsamples),[],[],10);

CELLS.matrixcount=matCounts(goodsamples,:);
CELLS.geneames=ALLEXPRESSION.genename;
CELLS.name=c(goodsamples)';
CELLS.pos=[ALLEXPRESSION.loc(goodsamples,2),ALLEXPRESSION.loc(goodsamples,1)];
CELLS.image=ALLEXPRESSION.image(sam);
CELLS.scale=ALLEXPRESSION.scale(sam);
%ALLCELLS{sam}=CELLS;
ALLCELLS = setfield(ALLCELLS,['sample',num2str(sam)],CELLS);
end
end