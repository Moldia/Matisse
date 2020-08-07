function[] = ISS_heatmap_cluster(CELLS)

CA=[];
DE=unique(CELLS.Clustering_ID);
for cant=1:size(unique(gro),1)
disp(cant)
group=DE(cant);
exp=mean(CELLS.expressionmatrix(ismember(gro,group),:),1);
CA=[CA;exp];
end
CAN=CA./mean(CA);
TO=CELLS.expressionnames;

figure;
heatmap(CAN,'Colormap',hot,'XDisplayLabels',TO,'YDisplayLabels',DE);

end