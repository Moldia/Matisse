function[RESULT] = heatmap_correlation_INPUTPRED(CELLS,MeanClassExp,ClassNames,GeneNames)

CA=[];
DE=unique(CELLS.Clustering_ID);
for cant=1:size(DE,1)
disp(cant)
group=DE(cant);
exp=mean(CELLS.expressionmatrix(ismember(CELLS.Clustering_ID,group),:),1);
CA=[CA;exp];
end
CAN=CA./mean(CA);
TO=CELLS.expressionnames;

figure;
heatmap(CAN,'Colormap',CELLS.colormap,'XDisplayLabels',TO,'YDisplayLabels',DE,'GridVisible','off',...
    'ColorScaling','scaledcolumns');
figure;

RESULT=zeros(size(CAN,1),size(MeanClassExp,1));
for col=1:size(CAN,1)
for row=1:size(MeanClassExp,1)
  vect=corrcoef(MeanClassExp(row,:),CAN(col,:));
  RESULT(col,row)=vect(1,2);
end

end

figure
heatmap(RESULT,'Colormap',CELLS.colormap,'XDisplayLabels',ClassNames,'YDisplayLabels',DE,'GridVisible','off')
figure
clustergram(RESULT(:,1:end-1),'Colormap',CELLS.colormap,'Rowlabels',DE,'Columnlabels',ClassNames(1:end-1),'Standardize','row')

end