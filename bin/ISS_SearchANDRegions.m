% search regions where all specified genes co-occur
% Xiaoyan, 2017
function []=ISS_SearchANDRegions(SPOTS,specific_genes)

%%
pos=SPOTS.pos;
name=SPOTS.name;
image=SPOTS.image;
scale=SPOTS.scale;

[name, pos] = removereads(name, 'NNNN', pos);
figure;
plotall(name, pos, image, scale)


% cooccur = search_reads_cooccur(name, pos, 100);
 cooccur = search_reads_cooccur(name, pos, 100,specific_genes);
% cooccur = search_reads_combs(name, pos, 500, 3, specific_genes);

end
