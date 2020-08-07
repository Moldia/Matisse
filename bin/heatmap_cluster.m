%Function: heatmap_cluster
%2020, Nilsson Lab, Scilifelab
%Objective: Show the expression profile of each of the clusters defined
%Input: Matisse object (BINS/CELLS)
%Output: Heatmap showing the expression of each cluster both clustegram and
%normal heatmap
function[] = heatmap_cluster(CELLS)
CA=[];
DE=unique(CELLS.Clustering_ID)';
for cant=1:size(DE,1)
disp(cant)
group=DE(cant);
exp=mean(CELLS.expressionmatrix(ismember(CELLS.Clustering_ID,group),:),1);
CA=[CA;exp];
end
CAN=CA./mean(CA);
TO=CELLS.expressionnames;

figure;
heatmap(CAN,'Colormap',CELLS.colormap,'XDisplayLabels',TO,'YDisplayLabels',DE,'GridVisible','off'); % 'ColorScaling','scaledcolumns'
figure;
clustergram(CAN,'Columnlabels',TO,'Rowlabels',DE,...
    'Colormap',CELLS.colormap)

end