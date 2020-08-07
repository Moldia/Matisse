function [SINGLE_CELL]=load_scRNAseq(path);
SINGLE_CELL=matisseMOD;
SC=readtable(path);
SINGLE_CELL.expressionmatrix=table2array(SC(:,2:end))';
SINGLE_CELL.expressionnames=table2cell(SC(:,1));
SINGLE_CELL.spotname=SC.Properties.VariableNames(2:end)';

end