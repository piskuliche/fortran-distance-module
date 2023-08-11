
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


def water_bond_list(nwaters):
    output = []
    ow = 1
    for water in range(nwaters):
        ow = 3*water + 1
        hw1 = ow + 1
        hw2 = ow + 2
        output.append([ow, hw1])
        output.append([ow, hw2])
    return np.array(output, dtype=int)


def water_osc_groups(nwaters):
    output = []
    count = 1
    for water in range(nwaters):
        output.append(count)
        output.append(count)
        output.append(count)
        count += 1
    return np.array(output, dtype=int)


if __name__ == "__main__":

    import MDAnalysis as mda
    import numpy as np

    u = mda.Universe("step5_1.tpr")

    aa = u.select_atoms("all")

    charges = aa.charges

    np.savetxt("charges.dat", charges, fmt="%10.5f")
    nwaters = int(len(charges)/3)

    oscs = water_hatom_locs(nwaters)
    bonds = water_bond_list(nwaters)
    osc_grps = water_osc_groups(nwaters)
    np.savetxt("oscs.dat", oscs, fmt="%d")
    with open("bonds.dat", 'w') as f:
        f.write("%d \n" % len(bonds))
        for bond in bonds:
            f.write("%d %d \n" % (bond[0], bond[1]))

    np.savetxt("osc_groups.dat", osc_grps, fmt="%d")
