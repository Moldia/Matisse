%Function: ROI_draw_onImage
%2020, Nilsson Lab, Scilifelab
%Objective: Select a region of interest for further analysis
%Input: Matisse object + sample selected (eg. 1)  
%    Optional input: for draw based on density, add a gene name on a list.
%Output: Matisse object including only the ROI selected 

function [SELECTIONS] = Draw_ROI(SPOTS,varargin)
if nargin==1
imagesel=1;
else
 imagesel=varargin{1}
end

roi_number=1;
allsamples=unique(SPOTS.sample);
SELECTED=ismember(SPOTS.sample,allsamples(imagesel));
SELECTIONS=SPOTS;
SELECTIONS.spotname=SELECTIONS.spotname(SELECTED,:);
SELECTIONS.location=SELECTIONS.location(SELECTED,:);
SELECTIONS.sample=SELECTIONS.sample(SELECTED,:);
SELECTIONS.images=SELECTIONS.images(imagesel);
SELECTIONS.scale=SELECTIONS.scale(imagesel);

name=SELECTIONS.spotname;
pos=SELECTIONS.location;
background_image=SELECTIONS.images{1};
scale=SELECTIONS.scale{1};



drawnow;
col = {'white' 'yellow' 'blue' 'tomato' 'palevioletred' 'bisque' 'pink' 'mediumspringgreen'};



[uNames, ~, iName] = unique(name);
cName = hist(iName, 1:numel(uNames));


disp('loading image..')
image = imread(background_image);

posScaled = correctcoord(pos, scale);

if nargin<3
  [Coord, Coord_write, blank] = drawpolygon(image, roi_number, scale);  
else
  if ismember(varargin(2),'spots')
  gene_density=varargin(2);
  bandwidth=SPOTS.ROI_bandwidth;
  density = gene_kde(name, pos, gene_density{1}, bandwidth, background_image, scale);
  [Coord, Coord_write, blank] = drawpolygon_spots(SELECTIONS, roi_number);    
  else
  gene_density=varargin(2);
  bandwidth=SPOTS.ROI_bandwidth;
  density = gene_kde(name, pos, gene_density{1}, bandwidth, background_image, scale);
  [Coord, Coord_write, blank] = drawpolygon_density(density, roi_number);
  end
  Coord_write(:,2:3) = Coord_write(:,2:3)./scale;
  coordinX=Coord{2,1};
  Coord{2,1}=coordinX./scale;
  coordinY=Coord{3,1};
  Coord{3,1}=coordinY./scale;
end
% create ROIs and get polygon coordinates (in original scale)








% remove empty polygons
if ~isempty(blank)
    Coord(:,blank) = [];
end



%%SAVE ROI where they are included
ROI_exist = 1:roi_number;
for se=1:roi_number 
    ROI_exist(blank) = [];
    [ROI_count, ROI_freq, ROI_proportion, ROI_area, In] = ...
    count_in_polygon(1,Coord(:,se),...
    pos, uNames, iName, cName);
    INROI(:,se)=In;
end


SELECTIONS.spotname=SELECTIONS.spotname(In,:);
SELECTIONS.location=SELECTIONS.location(In,:);
SELECTIONS.sample=SELECTIONS.sample(In,:);
if size(SELECTIONS.expressionmatrix,1)>1
SELECTIONS.expressionmatrix=SELECTIONS.expressionmatrix(In,:);
end
nameImag=SELECTIONS.images{1};
Imag=imread(nameImag);
Imag=imresize(Imag,1/SELECTIONS.scale{1});
Imag=Imag(min(Coord{3,1}):max(Coord{3,1}),min(Coord{2,1}):max(Coord{2,1}),:);
prompt = 'Inset ROI appendix for file:';
append = input(prompt)
newImName=[nameImag(1:end-4),'_',append,'.tif'];
imwrite(Imag,newImName);
SELECTIONS.images{1}=newImName;
SELECTIONS.scale{1}=1;
disp('A new image for the ROI will be created in the location of your original tif');
SELECTIONS.location(:,1)=SELECTIONS.location(:,1)-min(Coord{2,1});
SELECTIONS.location(:,2)=SELECTIONS.location(:,2)-min(Coord{3,1});

figure
imshow(imread(newImName)*20)
hold on
gscatter(SELECTIONS.location(:,1),SELECTIONS.location(:,2),SELECTIONS.spotname)
hold off

% count transcripts within polygons
ROI_exist = 1:roi_number;
ROI_exist(blank) = [];
[ROI_count, ROI_freq, ROI_proportion, ROI_area, In] = ...
    count_in_polygon(length(ROI_exist), Coord,...
    pos, uNames, iName, cName);

% % plot reads
% %plotall(name(In), pos(In,:), background_image, scale);
% %drawnow;
% 
% % plot polygons
% if length(col)<roi_number
%     col = repmat(col, 1, ceil(roi_number/length(col)));
% end
% figure;
% imshow(image,[]);
% hold on;
% for i = 1:length(ROI_exist)
%     plot(Coord{2,i}*scale, Coord{3,i}*scale,...
%         'LineWidth', 2, 'Color', rgb(col{ROI_exist(i)}));
% end
% legend(Coord(1,:), 'Location', 'NorthEastOutside', 'Color', [.6 .6 .6]);
% axis off;
% drawnow;
% 
% % bar plot
% figure;
% bh = bar(ROI_freq);
% set(gca, 'XTick', 1:numel(uNames), 'XTickLabel', uNames,...
%     'XLim', [0 numel(uNames)+1], 'XTickLabelRotation', 90,...
%     'FontSize', 6);
% legend(Coord(1,:), 'Location', 'NorthEastOutside', 'color', [.6 .6 .6]);
% ylabel('relative frequency');
% box off;
% for i = 1:length(bh)
%     set(bh(i), 'FaceColor', rgb(col{ROI_exist(i)}), 'EdgeColor',[.2 .2 .2], 'LineWidth',0.1);
% end
% drawnow;

% % write output file
% write_roi_countfile(fullfile(output_directory, 'ImageROICounts.csv'),...
%     Coord, uNames, ROI_count, ROI_area)

end
