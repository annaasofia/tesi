import numpy as np

SCALE = np.array([1e-3, 1e-3, 1e-6, 1e-6])

# order: [x, y, theta_x, theta_y]
COV_LIST = [
    np.array([  # run 1
        [ 3.82592775e+00,  7.08532402e-02,  2.94087313e+01,  2.44291188e+00],
        [ 7.08532402e-02,  7.22218856e+00, -6.30959887e+00,  2.81295028e+01],
        [ 2.94087313e+01, -6.30959887e+00,  8.68147418e+02,  6.94685664e+01],
        [ 2.44291188e+00,  2.81295028e+01,  6.94685664e+01,  1.91761445e+03],

    ]),
    np.array([  # run 2
        [ 3.92930430e+00,  9.40747550e-02,  3.02849998e+01,  2.86630787e+00],
        [ 9.40747550e-02,  7.23729626e+00, -6.11794240e+00,  2.81915332e+01],
        [ 3.02849998e+01, -6.11794240e+00,  8.73052463e+02,  7.20042455e+01],
        [ 2.86630787e+00,  2.81915332e+01,  7.20042455e+01,  1.91393472e+03],
    ]),
    np.array([  # run 3
        [ 4.27361685e+00, -8.57314662e-02,  3.73234726e+01,  1.69094804e+00],
        [-8.57314662e-02,  8.56729238e+00, -7.19990325e+00,  4.67594075e+01],
        [ 3.73234726e+01, -7.19990325e+00,  1.02760431e+03,  3.62442702e+01],
        [ 1.69094804e+00,  4.67594075e+01,  3.62442702e+01,  2.25058751e+03],
    ]),
    np.array([  # run 4
        [ 4.11661836e+00, -7.12587559e-02,  3.61415571e+01,  1.94711500e+00],
        [-7.12587559e-02,  8.09231439e+00, -6.69112700e+00,  4.10515883e+01],
        [ 3.61415571e+01, -6.69112700e+00,  1.01637172e+03,  4.18951958e+01],
        [ 1.94711500e+00,  4.10515883e+01,  4.18951958e+01,  2.18277297e+03],
    ]),
    np.array([  # run 5
        [ 4.15081002e+00, -6.26525519e-02,  3.63725329e+01,  2.00733617e+00],
        [-6.26525519e-02,  8.07673559e+00, -6.59538011e+00,  4.09247480e+01],
        [ 3.63725329e+01, -6.59538011e+00,  1.01900218e+03,  4.20490001e+01],
        [ 2.00733617e+00,  4.09247480e+01,  4.20490001e+01,  2.18161661e+03],
    ]),
]

# [mu_x, mu_y, mu_thx, mu_thy] in mm, mm, urad, urad
MU_LIST = [
    np.array([-0.21,  0.60, -0.93,  6.56]),   # run 1
    np.array([-0.35,  0.66, -1.60,  7.95]),   # run 2
    np.array([-0.60,  0.55, -2.79,  5.64]),   # run 3
    np.array([-0.42, -0.02, -1.73, -1.51]),   # run 4
    np.array([-0.38, -0.04, -1.37, -1.24]),   # run 5
]

N_LIST = np.array([12228420, 11094847, 8930171, 9437896, 10737344])


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


def cov_matrix(sigma_det_x=0.0, sigma_det_y=0.0):
    _, cov_avg = combine()
    cov_avg = cov_avg * np.outer(SCALE, SCALE)

    cov_avg[2, 2] -= sigma_det_x**2
    cov_avg[3, 3] -= sigma_det_y**2

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