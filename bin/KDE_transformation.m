% plotting on top of the density estimate
% Xiaoyan, 2020
function SPOTS = KDE_transformation(SPOTS)

%% modify here
sam=unique(SPOTS.sample);

for eachsam=1:size(sam,1)

sampl=sam(eachsam);
SELECTED=ismember(SPOTS.sample,sampl);
name=SPOTS.spotname(SELECTED,:);
pos=SPOTS.location(SELECTED,:);
image=SPOTS.images{eachsam};
scale=SPOTS.scale{eachsam};
%gene_density=SPOTS.Genes_of_interest{1};
bandwid=SPOTS.bandwidth;
% kde

SPO=unique(SPOTS.spotname);


AL=[];

for element=1:size(SPO,1)
disp(element);
gene_density=SPO(element);
bandwid=SPOTS.bandwidth;
density = gene_kde(name, pos, gene_density, bandwid, image, scale);
ALE=density;
V = ALE(:);
[X,Y]=find(ALE);
if size(AL,1)==0
    AL=[X,Y,V];
else
    AL=[AL,V]; 
end
end


[UMAPOUT,n,clusters]=run_umap(AL(1:30:end,3:end));









end






end 
