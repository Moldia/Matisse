%Function: OverlayTwoDensities
%2020, Nilsson Lab, Scilifelab
%Objective: Represent the KDE of two different genes
%Input: Matisse object (Need to have the genes to represents in
%Matisse.GenesofInterest
%Output: Image showing the KDE of two genes

function[SPOTS] = OverlayTwoDensities(SPOTS);
drawnow;
%% modify here
sam=unique(SPOTS.sample);

for eachsam=1:size(sam,1)
sampl=sam(eachsam);
SELECTED=ismember(SPOTS.sample,sampl);
name=SPOTS.spotname(SELECTED,:);
pos=SPOTS.location(SELECTED,:);
image=SPOTS.images{eachsam};
scale=SPOTS.scale{eachsam};
genes_density=SPOTS.OverlayTwoDensities_Genes_of_interest(:);
bandwidth=SPOTS.OverlayTwoDensities_bandwidth;


% unique transcripts
[uNames, ~, iName] = unique(name);
[p, q] = hist(iName, 1:length(uNames));

% image size
img = imread(image);
imsize = [size(img,1), size(img,2)];

layer = [1 0 0; 0 1 0;0 0 1];

% densities
I = zeros(ceil(imsize(1)/5), ceil(imsize(2)/5), 3, 'double');
for i = 1:2
    density = gene_kde(name, pos, genes_density{i}, bandwidth, image, scale);
    density = (density - min(density(:)))/max(density(:));   
    I = I + cat(3, density*layer(i,1), density*layer(i,2), density*layer(i,3));
end

figure;
imshow(I/max(I(:))*2);
title([genes_density{1} ' - ' genes_density{2} ' densities-' sampl{:} ]);

end

end