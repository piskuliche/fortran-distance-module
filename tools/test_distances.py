import numpy as np
import MDAnalysis as mda


L = [12.41380, 12.41380, 12.41380]
u = mda.Universe("../test_files/test.xyz")
atoms = u.select_atoms("all")

rc = 3
rc_sq = rc**2.
pos = atoms.positions
natoms = len(pos)

count = 0
for i in range(natoms):
    for j in range(i+1, natoms):
        dr = pos[i]-pos[j]
        dr = dr - L*np.round(dr/L)
        drsq = np.dot(dr, dr)
        if drsq < rc_sq:
            count += 1

print(np.max(pos))
print(np.min(pos))

print(count)
