function [varargout]=PlotTotalCounts(SPOTS,varargin);

[TABLE,ARG,AL,lab]=crosstab(SPOTS.spotname,SPOTS.sample);
figure
s=heatmap(lab(1:size(unique(SPOTS.sample),1),2),lab(1:size(unique(SPOTS.spotname),1),1),TABLE./sum(TABLE),'GridVisible','off','Colormap',SPOTS.colormap...
   , 'ColorScaling','scaledrows');
s.YLabel=('Genes');
s.XLabel=('Samples');
title('Total number of reads found on each sample');
figure
h=bar(categorical(lab(1:size(unique(SPOTS.spotname),1),1)),TABLE);
ylabel('Total counts');
xlabel('Genes');
set(h, {'DisplayName'}, lab(1:size(unique(SPOTS.sample),1),2))
title('Total number of reads found on each sample');


figure
h=bar(categorical(lab(1:size(unique(SPOTS.spotname),1),1)),TABLE./sum(TABLE));
set(h, {'DisplayName'}, lab(1:size(unique(SPOTS.sample),1),2))
ylabel('Percentage of reads');
xlabel('Genes');
title('Perceentage of genes on each sample');

tablexp=array2table(TABLE);
tablexp.Properties.VariableNames=lab(1:size(unique(SPOTS.sample),1),2);
tablexp.Properties.RowNames=lab(1:size(unique(SPOTS.spotname),1),1);

varargout{1}=tablexp;

end