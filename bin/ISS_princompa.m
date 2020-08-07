function[NEWMAT]= ISS_princompa(EXPRESSIONMAT)
[coeff,score,latent,~,explained] = pca(EXPRESSIONMAT.exp);
bar(explained)
title('Explained variability by different PC')
prompt = 'Number of PCA';
x = input(prompt)
NEWMAT=EXPRESSIONMAT;
NEWMAT.exp=score(:,1:x)*coeff(1:x,:);

end