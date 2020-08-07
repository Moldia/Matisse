function [MEANDIST] =densities_correlation(CELLSHIERARCHICAL,sam)
CLASSES=unique(CELLSHIERARCHICAL.spotname);
MEANDIST=zeros(size(CLASSES,1),size(CLASSES,1));
for START=1:size(CLASSES,1)
   disp(START);
    STARTCLASS=CLASSES(START);
   SUBSET.spotname=CELLSHIERARCHICAL.spotname(ismember(CELLSHIERARCHICAL.spotname,STARTCLASS),:);
   SUBSET.location=CELLSHIERARCHICAL.location(ismember(CELLSHIERARCHICAL.spotname,STARTCLASS),:); 
   for FINISH=1:size(CLASSES,1)
       FINISHCLASS=CLASSES(FINISH);
       SUBSET2.spotname=CELLSHIERARCHICAL.spotname(ismember(CELLSHIERARCHICAL.spotname,FINISHCLASS),:);
       SUBSET2.location=CELLSHIERARCHICAL.location(ismember(CELLSHIERARCHICAL.spotname,FINISHCLASS),:); 
       density = gene_kde(SUBSET.spotname, SUBSET.location,STARTCLASS, CELLSHIERARCHICAL.bandwidth, CELLSHIERARCHICAL.images{sam}, CELLSHIERARCHICAL.scale{sam});
       density2 = gene_kde(SUBSET2.spotname, SUBSET2.location,FINISHCLASS, CELLSHIERARCHICAL.bandwidth, CELLSHIERARCHICAL.images{sam}, CELLSHIERARCHICAL.scale{sam});
       MEANDIST(START,FINISH)=corr2(density,density2); 
       %DISTANCES=[];
       %for i=1:size(SUBSET.location,1)
       %   DISTANCES=[DISTANCES,mean((SUBSET.location(i,1)-SUBSET2.location(:,1)).^2 +(SUBSET.location(i,2)-SUBSET2.location(:,2)).^2)];
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
