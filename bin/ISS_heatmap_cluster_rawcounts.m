function[] = ISS_heatmap_cluster_rawcounts(EXPRESSIONMAT,gro)

CA=[];
DE=unique(gro);
for cant=1:size(unique(gro),1)
disp(cant)
group=DE(cant);
exp=mean(EXPRESSIONMAT.exp(ismember(gro,group),:),1);
CA=[CA;exp];
end
CAN=CA;
TO=EXPRESSIONMAT.genename;

figure;
heatmap(CAN,'Colormap',hot,'XDisplayLabels',TO,'YDisplayLabels',DE);

end