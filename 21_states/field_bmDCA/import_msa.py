import numpy as np

def import_msa_bmDCA(path_file, L):
    J2 = np.zeros((L, L, 21, 21), dtype=np.float64)
    h = np.zeros((L, 21), dtype=np.float64)
    with open(path_file, "r") as f:
        for line in f:
            l = line.rstrip("\n").split(" ")
            val = float(l[-1]) 
            if l[0] == "J":
                J2[int(l[1]), int(l[2]), int(l[3]), int(l[4])] = val
            elif l[0] == "h":
                h[int(l[1]), int(l[2])] = val
    # Symmetrize J2
    for i in range(J2.shape[0]):
        for j in range(i):
            J2[i, j, ...] = J2[j, i, ...].T
    return h,J2


def import_msa_bmDCA_dynamic(path_file):
    J2_idxs = []
    J2_vals = []
    h_idxs = []
    h_vals = []
    
    with open(path_file, "r") as f:
        for line in f:
            l = line.rstrip("\n").split(" ")
            val = float(l[-1])
            if l[0] == "J":
                J2_idxs.append([int(l[1]), int(l[2]), int(l[3]), int(l[4])])
                J2_vals.append(val)
            elif l[0] == "h":
                h_idxs.append([int(l[1]), int(l[2])])
                h_vals.append(val)

    J2_idxs = np.array(J2_idxs)
    h_idxs = np.array(h_idxs)
    
    L = np.max(h_idxs[:, 0]) + 1  # Infer L dynamically
    J2 = np.zeros((L, L, 21, 21), dtype=np.float64)
    h = np.zeros((L, 21), dtype=np.float64)

    h[tuple(h_idxs.T)] = h_vals
    J2[tuple(J2_idxs.T)] = J2_vals

    # Symmetrize J2
    for i in range(J2.shape[0]):
        for j in range(i):
            J2[i, j, ...] = J2[j, i, ...].T

    return h, J2
