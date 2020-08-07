function [] =ISS_heterogeneity_map(EXPRESSIONMAT,SPOTS)
TEST=[]; 
for s=1:size(EXPRESSIONMAT.exp,1)
TEST=[TEST,std(EXPRESSIONMAT.exp(s,:))];     
end
figure
%imshow(imresize(imread(SPOTS.image),1));
%hold on

scatter(EXPRESSIONMAT.loc(:,2),EXPRESSIONMAT.loc(:,1),10,sqrt(sqrt(sqrt(TEST))),'filled');

end
