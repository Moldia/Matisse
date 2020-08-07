function[NEWMAT,correlat]= exploreDEgenes(CELLS,varargin)

if nargin==1
imagesel=1;
else
 imagesel=varargin{1}
end


allsamples=unique(CELLS.sample);
SELECTED=ismember(CELLS.sample,CELLS.sample);

SELECTIONS=CELLS;
SELECTIONS.spotname=SELECTIONS.spotname(SELECTED,:);
SELECTIONS.location=SELECTIONS.location(SELECTED,:);
SELECTIONS.sample=SELECTIONS.sample(SELECTED,:);
SELECTIONS.images=SELECTIONS.images;
SELECTIONS.scale=SELECTIONS.scale(imagesel);
SELECTIONS.expressionmatrix=SELECTIONS.expressionmatrix(SELECTED,:);

GENES=SELECTIONS.expressionnames;

NM=[];
SM=[];
SD=[];
MN=[];
for g=1:size(GENES)
GE=GENES(g);
NM=[NM,GE];
SM=[SM,sum(SELECTIONS.expressionmatrix(:,g))];
SD=[SD,std(SELECTIONS.expressionmatrix(:,g))];
MN=[MN,mean(SELECTIONS.expressionmatrix(:,g))]
end

figure
scatter(MN,SD,10,'green','filled');
hold on
text(MN,SD,NM);
% set(gca,'xscale','log')
% set(gca,'yscale','log')
xlabel('Mean expression of each gene within bins/cells (logscale)')
ylabel('Standard deviation of each gene within bins/cells (logscale)')
title ('Number of spots VS standard deviation');


figure
bar(categorical(NM),SD./MN);
xlabel('Genes')
ylabel('Std/Mean')
end