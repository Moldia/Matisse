% draw polygons to create ROI 
% Polygons can be adjusted after they are created. When everything is done,
% double click to confirm the polygons
% Xiaoyan, 2018
function [INROI,Coord_write] = ISS_ROI_draw_onImage(SPOTS,output_directory,roi_number,varargin)


name=SPOTS.name;
pos=SPOTS.pos;
background_image=SPOTS.image;
scale=SPOTS.scale;



drawnow;
col = {'white' 'yellow' 'blue' 'tomato' 'palevioletred' 'bisque' 'pink' 'mediumspringgreen'};



[uNames, ~, iName] = unique(name);
cName = hist(iName, 1:numel(uNames));


disp('loading image..')
image = imread(background_image);

posScaled = correctcoord(pos, scale);

if nargin<4
  [Coord, Coord_write, blank] = drawpolygon(image, roi_number, scale);  
else
  gene_density=varargin(1);
  bandwidth=cell2mat(varargin(2));
  density = gene_kde(name, pos, gene_density{1}, bandwidth, background_image, scale);
  [Coord, Coord_write, blank] = drawpolygon_density(density, roi_number);
  Coord_write(:,2:3) = Coord_write(:,2:3)/scale;
end
% create ROIs and get polygon coordinates (in original scale)








% remove empty polygons
if ~isempty(blank)
    Coord(:,blank) = [];
end

% write coordinates to file
mkdir(output_directory);
fid = fopen(fullfile(output_directory, 'ImageROICoordinates.csv'), 'w');
fprintf(fid,'Polygon id,x coordiates,y coordinates\n');
fprintf(fid,'%d,%d,%d\n', Coord_write');
fclose(fid);

%%SAVE ROI where they are included
ROI_exist = 1:roi_number;
for se=1:roi_number 
    ROI_exist(blank) = [];
    [ROI_count, ROI_freq, ROI_proportion, ROI_area, In] = ...
    count_in_polygon(1,Coord(:,se),...
    pos, uNames, iName, cName);
    INROI(:,se)=In;
end



% count transcripts within polygons
ROI_exist = 1:roi_number;
ROI_exist(blank) = [];
[ROI_count, ROI_freq, ROI_proportion, ROI_area, In] = ...
    count_in_polygon(length(ROI_exist), Coord,...
    pos, uNames, iName, cName);

% plot reads
plotall(name(In), pos(In,:), background_image, scale);
drawnow;

% plot polygons
if length(col)<roi_number
    col = repmat(col, 1, ceil(roi_number/length(col)));
end
figure;
imshow(image,[]);
hold on;
for i = 1:length(ROI_exist)
    plot(Coord{2,i}*scale, Coord{3,i}*scale,...
        'LineWidth', 2, 'Color', rgb(col{ROI_exist(i)}));
end
legend(Coord(1,:), 'Location', 'NorthEastOutside', 'Color', [.6 .6 .6]);
axis off;
drawnow;

% bar plot
figure;
bh = bar(ROI_freq);
set(gca, 'XTick', 1:numel(uNames), 'XTickLabel', uNames,...
    'XLim', [0 numel(uNames)+1], 'XTickLabelRotation', 90,...
    'FontSize', 6);
legend(Coord(1,:), 'Location', 'NorthEastOutside', 'color', [.6 .6 .6]);
ylabel('relative frequency');
box off;
for i = 1:length(bh)
    set(bh(i), 'FaceColor', rgb(col{ROI_exist(i)}), 'EdgeColor',[.2 .2 .2], 'LineWidth',0.1);
end
drawnow;

% write output file
write_roi_countfile(fullfile(output_directory, 'ImageROICounts.csv'),...
    Coord, uNames, ROI_count, ROI_area)

end
