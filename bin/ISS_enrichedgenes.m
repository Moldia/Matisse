function [GLOBAL]=ISS_enrichedgenes(CELLS,type1,type2)
  SELEC.name=CELLS.name(ismember(CELLS.name,type1));
  SELEC.pos=CELLS.pos(ismember(CELLS.name,type1),:);
  SELEC.matrixcount=CELLS.matrixcount(ismember(CELLS.name,type1),:);
  SELEC.geneames=CELLS.geneames;
  
  COMPARE.name=CELLS.name(ismember(CELLS.name,type2));
  COMPARE.pos=CELLS.pos(ismember(CELLS.name,type2),:);
  COMPARE.matrixcount=CELLS.matrixcount(ismember(CELLS.name,type2));
  COMPARE.geneames=CELLS.geneames;
  
  VEC=[];
  for gen=1:size(SELEC.name,1)
      VEC=[VEC,min(sqrt((SELEC.pos(gen,1)-COMPARE.pos(:,1)).^2 + (SELEC.pos(gen,2)-COMPARE.pos(:,2)).^2))];
  end
  SELEC.dist=VEC;
  CLOSE=SELEC.matrixcount(SELEC.dist<prctile(SELEC.dist,10),:);
  FAR=SELEC.matrixcount(SELEC.dist>prctile(SELEC.dist,0),:);
  
  VECI=[];
  for gg=1:size(CLOSE,2)
      VECI=[VECI,ranksum(CLOSE(:,gg),FAR(:,gg))];
  end
  
  
  ENRICHMENT=(mean(CLOSE)-mean(FAR))./((mean(FAR)));
  
  
  GLOBAL=table(CELLS.geneames,ENRICHMENT',VECI');
  
  
end