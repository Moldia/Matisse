function []=ISS_plotCellReads(SPOTS,CellMap)
figure
image(label2rgb(CellMap, 'jet', 'w', 'shuffle'));
hold on
gscatter(SPOTS.pos(:,1),SPOTS.pos(:,2),SPOTS.name(:,:));
end