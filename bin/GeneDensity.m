% plotting on top of the density estimate
% Xiaoyan, 2020
function SPOTS = GeneDensity(SPOTS)

%% modify here
sam=unique(SPOTS.sample);

for eachsam=1:size(sam,1)
figure;

for geneselected=1:size(SPOTS.GeneDensity_Genes_of_interest,2);
genesel=SPOTS.GeneDensity_Genes_of_interest{geneselected};
sampl=sam(eachsam);
SELECTED=ismember(SPOTS.sample,sampl);
name=SPOTS.spotname(SELECTED,:);
pos=SPOTS.location(SELECTED,:);
image=SPOTS.images{eachsam};
scale=SPOTS.scale{eachsam};
gene_density=genesel;
bandwid=SPOTS.GeneDensity_bandwidth;
% kde
density = gene_kde(name, pos, gene_density, bandwid, image, scale);

% plot
subplot(1,size(SPOTS.GeneDensity_Genes_of_interest,2),geneselected);
imshow(density,[]);
hold on;
colormap(gca, SPOTS.GeneDensity_colormap);
title([gene_density ' density']);
end
%link axis if plot is bigger than 1
if size(SPOTS.GeneDensity_Genes_of_interest,2)>1
linkaxes();
end

end

end 
