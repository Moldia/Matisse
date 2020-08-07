function [MEANDIST] =ISS_densities_correlation(CELLSHIERARCHICAL,bandwid)
CELLSHIERARCHICAL.name=cellstr(num2str(CELLSHIERARCHICAL.name'));
CLASSES=unique(CELLSHIERARCHICAL.name);
MEANDIST=zeros(size(CLASSES,1),size(CLASSES,1));
for START=1:size(CLASSES,1)
    STARTCLASS=CLASSES(START);
   SUBSET.name=CELLSHIERARCHICAL.name(ismember(CELLSHIERARCHICAL.name,STARTCLASS),:);
   SUBSET.pos=CELLSHIERARCHICAL.pos(ismember(CELLSHIERARCHICAL.name,STARTCLASS),:); 
   for FINISH=1:size(CLASSES,1)
       FINISHCLASS=CLASSES(FINISH);
       SUBSET2.name=CELLSHIERARCHICAL.name(ismember(CELLSHIERARCHICAL.name,FINISHCLASS),:);
       SUBSET2.pos=CELLSHIERARCHICAL.pos(ismember(CELLSHIERARCHICAL.name,FINISHCLASS),:); 
       density = gene_kde(SUBSET.name, SUBSET.pos,STARTCLASS, bandwid, CELLSHIERARCHICAL.image, CELLSHIERARCHICAL.scale);
       density2 = gene_kde(SUBSET2.name, SUBSET2.pos,FINISHCLASS, bandwid, CELLSHIERARCHICAL.image, CELLSHIERARCHICAL.scale);
       MEANDIST(START,FINISH)=corr2(density,density2); 
       %DISTANCES=[];
       %for i=1:size(SUBSET.pos,1)
       %   DISTANCES=[DISTANCES,mean((SUBSET.pos(i,1)-SUBSET2.pos(:,1)).^2 +(SUBSET.pos(i,2)-SUBSET2.pos(:,2)).^2)];
       %end
       %MEANDIST(START,FINISH)=prctile(DISTANCES,1);
   end
end
%MEANDIST=(log(MEANDIST));
%figure
%heatmap(MEANDIST,'Colormap',hot,'XDisplayLabels',CLASSES,'YDisplayLabels',CLASSES);
MEANDIST(isnan(MEANDIST))=0;
figure
cg=clustergram(MEANDIST, 'ColumnLabels', CLASSES,'RowLabels',CLASSES);
cg.Colormap = hot;
end
