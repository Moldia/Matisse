%Function: AssignDomain
%2020, Nilsson Lab, Scilifelab
%Objective: Assign a domain to each cell
%Input: Matisse object containing binned data (CELLS/ BINS)
%Output: Low dimensional image our tissue, representing each color a different
%expression pattern.

function[CELLS]=AssignDomain(CELLS,BINS)
%% do not modify
% import data
samplesCELLS=unique(CELLS.sample);
samplesBINS=unique(BINS.sample);
if size(samplesCELLS)~=size(samplesBINS)
  error('Number of samples in the CELLS and BINS object is different');
end
%Loop for each sample
CELLS.domain=zeros(size(CELLS.sample));
CELLS.domain=num2cell(CELLS.domain);
for eachsam=1:size(samplesCELLS,1)
SAMPLES=ismember(CELLS.sample,samplesCELLS(eachsam));
CELLSLOC=CELLS.location(SAMPLES,:);
SAMPLESBINS=ismember(BINS.sample,samplesBINS(eachsam));
BINSLOC=BINS.location(SAMPLESBINS,:);
BINSTAG=BINS.spotname(SAMPLESBINS,:);
DOMAINFOUND=[];
for cell=1:size(CELLSLOC,1)
  CELLANA=CELLSLOC(cell,:);
 VA=BINSTAG(find(((CELLANA(:,1)-BINSLOC(:,1)).^2 + (CELLANA(:,2)-BINSLOC(:,2)).^2)==...
     min((CELLANA(:,1)-BINSLOC(:,1)).^2 + (CELLANA(:,2)-BINSLOC(:,2)).^2)));
DOMAINFOUND=[DOMAINFOUND,VA];    
end
CELLS.domain(SAMPLES)=DOMAINFOUND;

figure
subplot(1,2,1);
imshow(imresize(imread(CELLS.images{eachsam}),1/CELLS.scale{eachsam})*5);
title('RCPs coloured depending on the domain');
hold on
gscatter(CELLSLOC(:,1),CELLSLOC(:,2),DOMAINFOUND');
subplot(1,2,2);
imshow(imresize(imread(CELLS.images{eachsam}),1/CELLS.scale{eachsam})*5);
title('RCPs coloured depending on the gene name');
hold on
gscatter(CELLSLOC(:,1),CELLSLOC(:,2),CELLS.spotname(SAMPLES));
linkaxes();

SP=CELLS.spotname(SAMPLES);
DF=DOMAINFOUND;
[TAB,X,y,LABELS]=crosstab(SP,DF);
%TAB=TAB;
TAB=TAB./sum(TAB);
figure
bar(categorical(LABELS(1:size(unique(DF),2),2)),TAB','stacked');
legend(LABELS(1:size(unique(SP),1),1))
xlabel('Domains');
ylabel('Percentage of each cell/gene');
title('Distribution of different cells/genes in domains');
end
end