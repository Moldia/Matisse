%Function: LowDimensionalRGB
%2020, Nilsson Lab, Scilifelab
%Objective: Represent in RGB way the expression profile of samples
%Input: Matisse object containing binned data (CELLS/ BINS)
%Output: Low dimensional image our tissue, representing each color a different
%expression pattern.

function [CELLS_FILTERED]=Filter_cells(CELLS,varargin)
%% do not modify
% import data
CELLS_FILTERED=CELLS;
disp(['Initial number of genes included:  ', num2str(size(CELLS.expressionmatrix,2))]);
disp(['Initial number of cells/bins included: ',num2str(size(CELLS.spotname,1))]);
disp('FILTERS APPLIED');
disp(['MIN COUNTS PER GENE = ',num2str(CELLS.FC_MIN_counts_gene)]);
disp(['MIN COUNTS PER CELL = ',num2str(CELLS.FC_MIN_counts_cell)]);
disp(['MAX COUNTS PER CELL = ',num2str(CELLS.FC_MAX_counts_cell)]);

SELECTED_GENE=sum(CELLS_FILTERED.expressionmatrix)>CELLS.FC_MIN_counts_gene;
CELLS_FILTERED.expressionnames=CELLS_FILTERED.expressionnames(SELECTED_GENE);
CELLS_FILTERED.expressionmatrix=CELLS_FILTERED.expressionmatrix(:,SELECTED_GENE);

disp('%%%%%%%%%% AFTER FILTERING %%%%%%%%%%%%');
disp(['Number of genes filtered out: ',num2str(sum(~SELECTED_GENE))]);

figure
hist(sum(CELLS_FILTERED.expressionmatrix'),400);
hold on
xline(CELLS.FC_MIN_counts_cell,'--r');
xline(CELLS.FC_MAX_counts_cell,'--g');
xlabel('Number of reads/cell')
ylabel('Amount of cells')

SELECTED_CELL=sum(CELLS_FILTERED.expressionmatrix')>CELLS.FC_MIN_counts_cell & ...
    sum(CELLS_FILTERED.expressionmatrix')<CELLS.FC_MAX_counts_cell;
CELLS_FILTERED.expressionmatrix=CELLS_FILTERED.expressionmatrix(SELECTED_CELL,:);
CELLS_FILTERED.sample=CELLS_FILTERED.sample(SELECTED_CELL);
CELLS_FILTERED.location=CELLS_FILTERED.location(SELECTED_CELL,:);
CELLS_FILTERED.spotname=CELLS_FILTERED.spotname(SELECTED_CELL);
disp(['Number of cells filtered out: ',num2str(sum(~SELECTED_CELL))]);

disp(['Final number of genes included:  ', num2str(size(CELLS_FILTERED.expressionmatrix,2))]);
disp(['Final number of cells/bins included: ',num2str(size(CELLS_FILTERED.spotname,1))]);

end