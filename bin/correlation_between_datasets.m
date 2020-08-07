%Function: correlation_between_datasets
%2020, Nilsson Lab, Scilifelab
%Objective: Correlate between two datasets
%Input: -Matisse object of CELLS (including CellMap)
%       -Matisse object of single cell expression matrix derived from sigle cell data modified to
%       match gSet format
%Output: Matisse object contatining the gene expression of each cell and
%the cell type matching each cell 

function [] = correlation_between_datasets(CELLS,SINGLE_CELL,varargin)

CA=[];
DE=unique(CELLS.spotname);
for cant=1:size(DE,1)
disp(cant)
group=DE(cant);
exp=mean(CELLS.expressionmatrix(ismember(CELLS.spotname,group),:),1);
CA=[CA;exp];
end
%CAN=CA./mean(CA);
TO=CELLS.expressionnames;

CELLS2=CELLS;
CELLS2.expressionmatrix=CA;
CELLS2.spotname=DE;

if size(CELLS2.expressionmatrix,2) ~= size(SINGLE_CELL.expressionmatrix,2)
CELSEL=(ismember(CELLS2.expressionnames,SINGLE_CELL.expressionnames));
CELLS2.expressionnames=CELLS2.expressionnames(CELSEL);
CELLS2.expressionmatrix=CELLS2.expressionmatrix(:,CELSEL);
[CELLS2.expressionnames,ORDER]=sort(CELLS2.expressionnames);
CELLS2.expressionmatrix=CELLS2.expressionmatrix(:,ORDER);

CELSEL2=(ismember(SINGLE_CELL.expressionnames,CELLS2.expressionnames));
SINGLE_CELL.expressionnames=SINGLE_CELL.expressionnames(CELSEL2);
SINGLE_CELL.expressionmatrix=SINGLE_CELL.expressionmatrix(:,CELSEL2);
[SINGLE_CELL.expressionnames,ORDER]=sort(SINGLE_CELL.expressionnames);
SINGLE_CELL.expressionmatrix=SINGLE_CELL.expressionmatrix(:,ORDER);
end

ALLTAB=array2table(SINGLE_CELL.expressionmatrix);
ALLTAB.Properties.VariableNames=SINGLE_CELL.expressionnames;
ALLTAB.Properties.RowNames=SINGLE_CELL.spotname;

EXPRESSEDCLUSTERS=array2table(CELLS2.expressionmatrix);
EXPRESSEDCLUSTERS.Properties.VariableNames=CELLS2.expressionnames;
EXPRESSEDCLUSTERS.Properties.RowNames=CELLS2.spotname;


corrmat=zeros(size(ALLTAB,1),size(EXPRESSEDCLUSTERS,1));
for x=1:size(ALLTAB,1)
   for y=1:size(EXPRESSEDCLUSTERS,1)
      coe=corrcoef(table2array(ALLTAB(x,:)),table2array(EXPRESSEDCLUSTERS(y,:)));
       corrmat(x,y)=coe(2,1);
   end
end


figure
clustergram(corrmat,'RowLabels',ALLTAB.Properties.RowNames,...
   'Colormap',CELLS.colormap,'ColumnLabels',EXPRESSEDCLUSTERS.Properties.RowNames,'Standardize','row');

figure
heatmap(EXPRESSEDCLUSTERS.Properties.RowNames,ALLTAB.Properties.RowNames,corrmat,'GridVisible','off','Colormap',parula,'ColorScaling',"scaledcolumns");

if nargin>2

EXP=table2array(EXPRESSEDCLUSTERS(ismember(EXPRESSEDCLUSTERS.Properties.RowNames,varargin{1}),:));
ALL=table2array(ALLTAB(ismember(ALLTAB.Properties.RowNames,varargin{2}),:));

figure
scatter(EXP,ALL,'.','red');
hold on 
text(EXP,ALL,EXPRESSEDCLUSTERS.Properties.VariableNames);
xlabel(['Number of spots (ISS)-',varargin{1}])
ylabel(['Number of spots (SCRNA)-',varargin{2}])
s=corrcoef(EXP,ALL);
title([' correlation between clusters== ',num2str(s(1,2))]);
end 


end
