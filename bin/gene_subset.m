%Function: gene_subset
%2020, Nilsson Lab, Scilifelab
%Objective: Select a subset of genes/cells
%samples into a single object
%Input: Matisse object
%    Optional input: cell list of genes to be selected
%Output: Matisse object including olny the gene selected 

function [SPOTS]=gene_subset(SPOTS,varargin)
%%%%%%If we are working with bins, just select the columns in expression
%%%%%%matrix
if size(SPOTS.expressionmatrix,1)>0;
if nargin<2
cNames=(unique(SPOTS.expressionnames));
cbValues = checkboxes(cNames);
cbValues=cell2table(cbValues);
cbValues=cbValues(cbValues.cbValues2==1,:);
GENES=cbValues.cbValues1;
else
GENES=varargin{1}; 
end
GENESELE=ismember(SPOTS.expressionnames,GENES);
SPOTS.expressionnames=SPOTS.expressionnames(GENESELE);
SPOTS.expressionmatrix=SPOTS.expressionmatrix(:,GENESELE);

%%%%% If we work with spots itself, select the spots 
else
if nargin<2
cNames=(unique(SPOTS.spotname));
cbValues = checkboxes(cNames);
cbValues=cell2table(cbValues);
cbValues=cbValues(cbValues.cbValues2==1,:);
GENES=cbValues.cbValues1;
else
   GENES=varargin{1}; 
end
GENESELE=ismember(SPOTS.spotname,GENES);
SPOTS.sample=SPOTS.sample(GENESELE,:);
SPOTS.spotname=SPOTS.spotname(GENESELE,:);
SPOTS.location=SPOTS.location(GENESELE,:);

end

end