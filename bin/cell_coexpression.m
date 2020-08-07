function cell_coexpression(CELLS,imagesel)

roi_number=1;
allsamples=unique(CELLS.sample);
SELECTED=ismember(CELLS.sample,allsamples(imagesel));
SELECTIONS=CELLS;
SELECTIONS.spotname=SELECTIONS.spotname(SELECTED,:);
SELECTIONS.location=SELECTIONS.location(SELECTED,:);
SELECTIONS.sample=SELECTIONS.sample(SELECTED,:);
SELECTIONS.images=SELECTIONS.images;
SELECTIONS.scale=SELECTIONS.scale(imagesel);
SELECTIONS.expressionmatrix=SELECTIONS.expressionmatrix(SELECTED,:)

COEXP=zeros(size(SELECTIONS.expressionnames,1),size(SELECTIONS.expressionnames,1));

for gene1=1:size(SELECTIONS.expressionnames);
 MATR=SELECTIONS.expressionmatrix(SELECTIONS.expressionmatrix(:,gene1)>0,:);
 for gene2=1:size(SELECTIONS.expressionnames);
    percent=sum(MATR(:,gene2)>0)/size(MATR,1); %Initially >0
    if gene1==gene2
        percent=1;
    end
    COEXP(gene1,gene2)=percent;
 end
end

figure
heatmap(SELECTIONS.expressionnames,SELECTIONS.expressionnames,COEXP,'Colormap',SELECTIONS.colormap,'Xlabel','Percentage of this gene','Ylabel','Gene analyzed',...
    'GridVisible','off','Title','Coexpression without normalization');

figure
heatmap(SELECTIONS.expressionnames,SELECTIONS.expressionnames,COEXP,'Colormap',SELECTIONS.colormap,'Xlabel','Percentage of this gene','Ylabel','Gene analyzed',...
    'ColorScaling','ScaledRows','GridVisible','off','Title','Coexpression with scaled rows');

end


