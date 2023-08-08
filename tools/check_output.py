def create_pairwise(a, b):
    pairs = []
    c0, c1, c2 = 0, 0, 0
    for i in range(len(a)):
        if a[i] == b[i]:
            c0 += 1
        elif a[i] > b[i]:
            c1 += 1
            pairs.append((min(a[i], b[i]), max(a[i], b[i])))
        else:
            c2 += 1
            pairs.append((min(a[i], b[i]), max(a[i], b[i])))
    print("a=b: ", c0)
    print("a>b: ", c1)
    print("a<b: ", c2)
    print("Total: ", c0 + c1 + c2)
    return pairs


def check_unique_pairs(a, b):
    pairs = create_pairwise(a, b)
    unique_pairs = set(tuple(pairs))
    print("Original pairs: ", len(pairs))
    print("Unique pairs: ", len(unique_pairs))
    return unique_pairs


def compare_output(a, b):
    c_anotinb = 0
    for pair in a:
        if pair not in b:
            # print(pair, "in a and not in b")
            c_anotinb += 1

    c_bnotina = 0
    for pair in b:
        if pair not in a:
            # print(pair, "in b and not in a")
            c_bnotina += 1

    print("a not in b: ", c_anotinb)
    print("b not in a: ", c_bnotina)

    return


def check_dists(filestr, a, b):
    drs = []
    try:
        drs = np.genfromtxt(filestr + "_dr.dat", dtype=float)
    except:
        raise FileNotFoundError("File dr is not a valid file")

    mloc = np.where(drs == max(drs))[0][0]
    print("max", max(drs), a[mloc], b[mloc])
    mloc = np.where(drs == min(drs))[0][0]
    print("min", min(drs), a[mloc], b[mloc])
    print("unique", len(np.unique(drs)))


if __name__ == "__main__":
    import sys
    import numpy as np

    # Parse the system arguments
    if len(sys.argv) != 3:
        print("Usage: python3 check_double_loop.py <filename>")
        raise ValueError("Invalid number of arguments")

    if type(sys.argv[1]) != str:
        print("Usage: python3 check_double_loop.py <filename>")
        raise ValueError("Filename 1 must be a string")

    if type(sys.argv[2]) != str:
        print("Usage: python3 check_double_loop.py <filename>")
        raise ValueError("Filename 2 must be a string")

    filestr1 = sys.argv[1]
    filestr2 = sys.argv[2]
    filestrs = [filestr1, filestr2]

    unique_pairs = {}
    for i in range(2):
        print("********************")
        print("Checking file: ", filestrs[i])
        # load the files into numpy arrays
        id1, id2 = [], []
        try:
            id1 = np.genfromtxt(filestrs[i]+"_id1.dat", dtype=int)
        except:
            raise FileNotFoundError("File 1 is not a valid file")

        try:
            id2 = np.genfromtxt(filestrs[i] + "_id2.dat", dtype=int)
        except:
            raise FileNotFoundError("File 2 is not a valid file")

        check_dists(filestrs[i], id1, id2)
        # Check the unique pairs in each file.
        tmp = check_unique_pairs(id1, id2)
        unique_pairs[i] = tmp
    print("********************")
    compare_output(unique_pairs[0], unique_pairs[1])
