% plotting on top of the density estimate
% Xiaoyan, 2020
function[] = ISS_GeneDensity(SPOTS, gene_density, bandwid, image, scale)

%% modify here
name=SPOTS.name;
pos=SPOTS.pos;
image=SPOTS.image;
scale=SPOTS.scale;

% kde
density = gene_kde(name, pos, gene_density, bandwid, image, scale);

% plot
figure;
imshow(density,[]);
hold on;
colormap(gca, parula);
title([gene_density ' density']);
end 
