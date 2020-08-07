function [GLOBAL]=enrichedgenes(CELLS,imagesel,type1,type2)
  roi_number=1;
allsamples=unique(CELLS.sample);
SELECTED=ismember(CELLS.sample,allsamples(imagesel));



SELECTIONS=CELLS;
SELECTIONS.spotname=SELECTIONS.spotname(SELECTED,:);
SELECTIONS.location=SELECTIONS.location(SELECTED,:);
SELECTIONS.sample=SELECTIONS.sample(SELECTED,:);
SELECTIONS.images=SELECTIONS.images{imagesel};
SELECTIONS.scale=SELECTIONS.scale(imagesel);


%cel




  SELEC.name=SELECTIONS.spotname(ismember(SELECTIONS.spotname,type1));
  SELEC.pos=SELECTIONS.location(ismember(SELECTIONS.spotname,type1),:);
  SELEC.matrixcount=SELECTIONS.expressionmatrix(ismember(SELECTIONS.spotname,type1),:);
  SELEC.geneames=SELECTIONS.expressionnames;
  
  COMPARE.name=SELECTIONS.spotname(ismember(SELECTIONS.spotname,type2));
  COMPARE.pos=SELECTIONS.location(ismember(SELECTIONS.spotname,type2),:);
  COMPARE.matrixcount=SELECTIONS.expressionmatrix(ismember(SELECTIONS.spotname,type2));
  COMPARE.geneames=SELECTIONS.expressionnames;
  
  VEC=[];
  for gen=1:size(SELEC.name,1)
      VEC=[VEC,min(sqrt((SELEC.pos(gen,1)-COMPARE.pos(:,1)).^2 + (SELEC.pos(gen,2)-COMPARE.pos(:,2)).^2))];
  end
  SELEC.dist=VEC;
  CLOSE=SELEC.matrixcount(SELEC.dist<prctile(SELEC.dist,10),:);
  FAR=SELEC.matrixcount(SELEC.dist>prctile(SELEC.dist,10),:);
  
  VECI=[];
  for gg=1:size(CLOSE,2)
      VECI=[VECI,ranksum(CLOSE(:,gg),FAR(:,gg))];
  end
  
  
  ENRICHMENT=((mean(CLOSE)-mean(FAR)))./((std([FAR;CLOSE])));
  
  
  GLOBAL=table(SELECTIONS.expressionnames,ENRICHMENT',VECI');
  GLOBAL.Properties.VariableNames={'Gene','Change in expression','pvalue'}
  
  figure
  gscatter(VECI',ENRICHMENT,SELECTIONS.expressionnames);
  hold on
  text(VECI',ENRICHMENT',SELECTIONS.expressionnames,'FontSize',8);
  xline(0.05,'-.','0.05 p-value','DisplayName','Average Sales');
  yline(0,'r-.','No diff. expr','DisplayName','Average Sales');
  stem(VECI',ENRICHMENT)
  ylabel('Enrichment (fold)');
  xlabel('p-value');
  title(['Genes enriched in cells from cluster ',type1,' close to cells from',type2]);
end