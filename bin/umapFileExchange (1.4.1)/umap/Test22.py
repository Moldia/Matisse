
import csv
import umap
import numpy as np
import matplotlib.pyplot as plt

with open('umap/sample55k.csv', 'r') as f:
    wines = list(csv.reader(f, delimiter=","))

inData = np.array(wines[1:], dtype=np.float)
inData2 = inData.copy()
inData2[22000,0] = inData2[22000,0] + 0.001

with open('umap/sampleBalbcLabeled55k.csv', 'r') as f:
    labels = list(csv.reader(f, delimiter=","))

reducer1 = umap.UMAP(verbose=True)
reducer1 = reducer1.fit(inData)
embedding1 = reducer1.embedding_
embedding2 = reducer1.transform(inData2)

labels = np.array(labels[1:], dtype=np.float)
labels = labels[:,10]


reducer2 = umap.UMAP(verbose=True)
reducer2 = reducer2.fit(inData, labels)
embedding3 = reducer2.embedding_
embedding4 = reducer2.transform(inData2)

plt.scatter(*embedding3.T, s=0.3, alpha=1.0)
plt.show()
plt.scatter(*embedding4.T, s=0.3, alpha=1.0)
plt.show()



