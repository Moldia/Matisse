%Function: LowDimensionalRGB
%2020, Nilsson Lab, Scilifelab
%Objective: Represent in RGB way the expression profile of samples
%Input: Matisse object containing binned data (CELLS/ BINS)
%Output: Low dimensional image our tissue, representing each color a different
%expression pattern.

function [CELLS_FILTERED]=Normalize_cells(CELLS,varargin)
%% do not modify
% import data
CELLS_FILTERED=CELLS;

if ismember(CELLS.NC_method,'log2')
CELLS_FILTERED.expressionmatrix=log2(CELLS_FILTERED.expressionmatrix+1);
end

if ismember(CELLS.NC_method,'log10')
CELLS_FILTERED.expressionmatrix=log10(CELLS_FILTERED.expressionmatrix);
end

if ismember(CELLS.NC_method,'standardize')
CELLS_FILTERED.expressionmatrix=(CELLS_FILTERED.expressionmatrix-mean(CELLS_FILTERED.expressionmatrix))./std(CELLS_FILTERED.expressionmatrix);
end


end