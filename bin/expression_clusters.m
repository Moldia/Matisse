%Function: gene_subset
%2020, Nilsson Lab, Scilifelab
%Objective: Select a subset of genes/cells
%samples into a single object
%Input: Matisse object
%    Optional input: cell list of genes to be selected
%Output: Matisse object including olny the gene selected 

function []=expression_clusters(SPOTS,method,varargin)
if nargin<3
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
SPOTS.expressionmatrix=SPOTS.expressionmatrix(GENESELE,:);

SI=unique(SPOTS.spotname);

if ismember(method,'violin') 
ES=repelem({'r','b','m','k','g','y','c',[0 0.5 1],[0 0 0.7],[0.1 0.7 0.7],[0.6 0 0],[0.4 0 0.3]},size(SPOTS.expressionmatrix,2));
 NES=randperm(numel(ES));
 ES=ES(NES);

for unn=1:8:size(SI,1)
S=SI(unn:min(unn+8,size(SI,1)));
figure
for Se=1:size(S,1)
 subplot(round(size(S,1)/2),2,Se);
 clust=S(Se);
 SELIS=ismember(SPOTS.spotname,clust);
 SEL=SPOTS.expressionmatrix(SELIS,:);
 title(clust);
 h= distributionPlot(SEL,'color',ES(1:size(SEL,2)),'xNames',SPOTS.expressionnames,'addSpread',0,'showMM',0);
 linkaxes();
 hold off
end

end

end
if ismember(method,'boxplot')
ES=repelem({'r','b','m','k','g','y','c',[0 0.5 1],[0 0 0.7],[0.1 0.7 0.7],[0.6 0 0],[0.4 0 0.3]},size(SPOTS.expressionmatrix,2));
 NES=randperm(numel(ES));
 ES=ES(NES);

for unn=1:8:size(SI,1)
S=SI(unn:min(unn+8,size(SI,1)));
figure
for Se=1:size(S,1)
 subplot(round(size(S,1)/2),2,Se);
 clust=S(Se);
 SELIS=ismember(SPOTS.spotname,clust);
 SEL=SPOTS.expressionmatrix(SELIS,:);
 title(clust);
 boxplot(SEL,SPOTS.expressionnames);
 linkaxes();
 hold off
end

end

end

end