
def water_hatom_locs(nwaters):
    """
    Return the indices of the hydrogen atoms in a water molecule.
    """
    output = []
    for water in range(nwaters):
        output.append(0)
        output.append(1)
        output.append(1)
    return np.array(output, dtype=int)


if __name__ == "__main__":

    import MDAnalysis as mda
    import numpy as np

    u = mda.Universe("step5_1.tpr")

    aa = u.select_atoms("all")

    charges = aa.charges

    np.savetxt("charges.dat", charges)
    nwaters = int(len(charges)/3)

    oscs = water_hatom_locs(nwaters)
    np.savetxt("oscs.dat", oscs, fmt="%d")
