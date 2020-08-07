% plotting on top of the density estimate
% Xiaoyan, 2018
function[] = ISS_PlotOnDensity(SPOTS, gene_density, bandwid)

%% modify here
name=SPOTS.name;
pos=SPOTS.pos;
image=SPOTS.image;
scale=SPOTS.scale;
% kde
density = gene_kde(name, pos, gene_density, bandwid, image, scale);

% plot
figure;
imshow(density, []);
colormap(gca, parula);
hold on;
plotall(name, correctcoord(pos, scale));
update_legend(gca, genes_to_plot);
title([gene_density ' density']);
