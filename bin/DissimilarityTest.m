% test the null hypothesis if two groups of data points are drawn from the same continuous distribution 
% two-sample Kolmogorov-Smirnov test
% Xiaoyan, 2017

function[] = DissimilarityTest(SPOTS,varargin);

% genes with read lower than this threshold will not be used in analysis


%SET THE VARARGIN VARIABLES
if (nargin>1)
    count_threshold=cell2mat(varargin(1));  
else
    count_threshold = SPOTS.count_threshold;  
end

if (nargin>2)
    bin_size=cell2mat(varargin(2));  
else
    bin_size = SPOTS.hexbin_size;  
end


sam=unique(SPOTS.sample);

for eachsam=1:size(sam,1)
sampl=sam(eachsam);
SELECTED=ismember(SPOTS.sample,sampl);
name=SPOTS.spotname(SELECTED,:);
pos=SPOTS.location(SELECTED,:);

[uNames, ~, idxName] = unique(name);
cNames = hist(idxName, 1:max(idxName));

% remove NNNN and any reads<count_threshold
[name, pos] = removereads(name, [uNames(cNames<count_threshold); {'NNNN'}], pos);
[uNames, ~, idxName] = unique(name);

pos_bin = ceil(pos/bin_size);
pos_bin = max(pos_bin(:,1))*(pos_bin(:,1)-1) + pos_bin(:,2);


ksstat = zeros(length(uNames), 2);
for i = 1:length(uNames)
    for j = i:length(uNames)
        [~, p, ks2stat] = kstest2(pos_bin(idxName==i), pos_bin(idxName==j),'tail','larger');
        ksstat(i,j,1) = ks2stat;
        ksstat(i,j,2) = p;
        ksstat(j,i,1) = ks2stat;
        ksstat(j,i,2) = p;
    end
end
% figure; imagesc(ksstat(:,:,2));
% set(gca, 'xtick', 1:length(uniName), 'xticklabel', uniName,...
%     'ytick', 1:length(uniName), 'yticklabel', uniName,...
%     'xticklabelrotation', 90);


% sorting
L = linkage(ksstat(:,:,2));
order = dendroperm(L);

figure; imagesc(ksstat(order,order,2)); grid on;
set(gca, 'xtick', 1:length(uNames), 'xticklabel', uNames(order),...
    'ytick', 1:length(uNames), 'yticklabel', uNames(order),...
    'xticklabelrotation', 90,'FontSize',10);

c = colorbar;
set(gca,'colormap',SPOTS.colormap);
c.Label.String = 'Kolmogorov-Smirnov test p-value';
c.Label.FontSize = 12;
%title([sampl{:} ]);

end
end
