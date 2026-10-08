import numpy as np

SCALE = np.array([1e-3, 1e-3, 1e-6, 1e-6])

# order: [x, y, theta_x, theta_y]
COV_LIST = [
    np.array([  # run 1
        [   4.74346824,    0.01606491,   89.41765855,    7.18923536],
        [   0.01606491,    8.72738918,   -9.9530266 ,   71.62957414],
        [  89.41765855,   -9.9530266 , 5362.42227158,  390.61366706],
        [   7.18923536,   71.62957414,  390.61366706, 4824.08887831],
    ]),
    np.array([  # run 2
        [   4.81758545,    0.0369429 ,   87.49673916,    7.70434782],
        [   0.0369429 ,    8.72286463,   -9.77122957,   69.79839735],
        [  87.49673916,   -9.77122957, 5109.67705187,  399.79783174],
        [   7.70434782,   69.79839735,  399.79783174, 4657.48262579],
    ]),
    np.array([  # run 3
        [   5.30335239,   -0.24781449,   96.73768344,    4.86882942],
        [  -0.24781449,   10.2478887 ,  -13.26363642,   92.85475119],
        [  96.73768344,  -13.26363642, 5139.62932306,  226.98102882],
        [   4.86882942,   92.85475119,  226.98102882, 5310.7596767 ],
    ]),
    np.array([  # run 4
        [   5.15235211,   -0.21145735,   95.7287431 ,    5.86400589],
        [  -0.21145735,    9.73438287,  -11.89446127,   85.77557724],
        [  95.7287431 ,  -11.89446127, 5141.51246685,  274.29155722],
        [   5.86400589,   85.77557724,  274.29155722, 5117.69360654],
    ]),
    np.array([  # run 5
        [   5.20988254,   -0.20068793,   98.20184138,    5.69484793],
        [  -0.20068793,    9.736861  ,  -11.81183585,   88.44745455],
        [  98.20184138,  -11.81183585, 5343.63042278,  248.74821146],
        [   5.69484793,   88.44745455,  248.74821146, 5357.36343426],
    ]),
]

# [mu_x, mu_y, mu_thx, mu_thy] in mm, mm, urad, urad
MU_LIST = [
    np.array([-0.21,  0.60, -0.93,  6.79]),   # run 1
    np.array([-0.35,  0.66, -1.61,  8.18]),   # run 2
    np.array([-0.60,  0.55, -2.79,  6.05]),   # run 3
    np.array([-0.42, -0.02, -1.74, -1.03]),   # run 4
    np.array([-0.38, -0.04, -1.38, -0.03]),   # run 5
]

N_LIST = np.array([12480454, 11322311, 9097215, 9617081, 10940060])


def combine():
    # cov_tot = sum_i w_i * ( C_i + (mu_i - mu_avg)(mu_i - mu_avg)^T )

    assert len(COV_LIST) == len(MU_LIST) == len(N_LIST)
    w = N_LIST / N_LIST.sum()

    mu_avg = sum(wi * m for wi, m in zip(w, MU_LIST))
    cov_avg = sum(
        wi * (C + np.outer(m - mu_avg, m - mu_avg))
        for wi, C, m in zip(w, COV_LIST, MU_LIST)
    )
    return mu_avg, cov_avg


def cov_matrix():
    _, cov_avg = combine()
    cov_avg = cov_avg * np.outer(SCALE, SCALE)
    eigvals = np.linalg.eigvalsh(cov_avg)
    return cov_avg, eigvals, N_LIST


def mean_vector():
    mu_avg, _ = combine()
    return mu_avg * SCALE


if __name__ == "__main__":
    mu, cov = combine()
    print("mu   [mm, mm, urad, urad]:", mu)
    print("sigma[mm, mm, urad, urad]:", np.sqrt(np.diag(cov)))
    print("min eigenvalue (SI):", cov_matrix()[1].min())