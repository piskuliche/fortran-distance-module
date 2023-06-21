import numpy as np

v1, a1, a2, nv1, na1, na2 = np.genfromtxt("map_reg.dat", usecols=(0,1,2,3,4,5), unpack=True)

i=0
for v in nv1:
    if v not in v1:
        print(v, na1[i], na2[i])
    i+=1

print(len(np.unique(v1)))
print(len(v1[v1!=0]))
print(len(np.unique(nv1)))
print(len(nv1[nv1!=0]))

import collections
print([item for item, count in collections.Counter(v1).items() if count>1])
print(np.where(v1==0.88673687))
print(v1[2], a1[2], a2[2])
print(v1[3], a1[3], a2[3])