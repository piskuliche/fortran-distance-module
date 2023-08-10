import MDAnalysis as mda


u = mda.Universe("step5_1.gro", "step5_1.xtc")
aa = u.select_atoms("all")

with mda.Writer("temp.xyz", aa.n_atoms) as W:
    f = open("L.dat", 'w')
    for ts in u.trajectory[:10]:
        W.write(aa)
        f.write("%10.5f %10.5f %10.5f\n" %
                (ts.dimensions[0], ts.dimensions[1], ts.dimensions[2]))
    f.close()
