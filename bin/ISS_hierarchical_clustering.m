function [CELLS]=ISS_hierarchical_clustering(EXPRESSIONMAT,SPOTS,n_clust)
matCounts=EXPRESSIONMAT.exp;
Z = linkage(matCounts(:,:),'ward');
c = cluster(Z,'Maxclust',n_clust);
%figure
%X=tsne(matCounts, 'NumDimensions', 3, 'NumPCAComponents', 4);
%scatter3(X(:,1),X(:,2),X(:,3),10,c)

im=imread(SPOTS.image);
ima=imresize(im,1/SPOTS.scale);
figure
imshow(ima*10);
hold on
g=gscatter(EXPRESSIONMAT.loc(:,2),EXPRESSIONMAT.loc(:,1),c,[],[],10);

CELLS.matrixcount=matCounts;
CELLS.geneames=EXPRESSIONMAT.genename;
CELLS.name=c;
CELLS.pos=[EXPRESSIONMAT.loc(:,2),EXPRESSIONMAT.loc(:,1)];
CELLS.image=SPOTS.image;
CELLS.scale=SPOTS.scale;
end