import numpy as np
from scipy.stats import norm


def cov_matrix():
    cov_run1 = np.array([
        [ 3.25259585e-01,  9.19410320e-03,  2.48120848e+00,  2.79216632e-01],
        [ 9.19410320e-03,  5.12463330e+00, -5.18479884e+00,  8.37145297e+00],
        [ 2.48120848e+00, -5.18479884e+00,  8.06634056e+02,  1.04718809e+02],
        [ 2.79216632e-01,  8.37145297e+00,  1.04718809e+02,  1.82409073e+03],
    ])
    n_run1 = 4263962

    cov_run2 = np.array([
        [ 3.25411796e-01,  9.97010963e-03,  2.49535577e+00,  2.98410762e-01],
        [ 9.97010963e-03,  5.16976789e+00, -5.28622937e+00,  8.68507428e+00],
        [ 2.49535577e+00, -5.28622937e+00,  8.04843024e+02,  1.05618050e+02],
        [ 2.98410762e-01,  8.68507428e+00,  1.05618050e+02,  1.82443400e+03],
    ])
    n_run2 = 3807353

    cov_run3 = np.array([
        [ 3.24311107e-01, -2.77618082e-03,  2.73816756e+00,  1.49477197e-01],
        [-2.77618082e-03,  5.99476224e+00, -4.98829806e+00,  2.10951692e+01],
        [ 2.73816756e+00, -4.98829806e+00,  8.56562329e+02,  7.47785904e+01],
        [ 1.49477197e-01,  2.10951692e+01,  7.47785904e+01,  2.12594105e+03],
    ])
    n_run3 = 2651593

    cov_run4 = np.array([
        [ 3.18082275e-01, -8.60134132e-03,  2.73927086e+00,  1.69471687e-01],
        [-8.60134132e-03,  5.84661715e+00, -5.67027042e+00,  1.99983894e+01],
        [ 2.73927086e+00, -5.67027042e+00,  9.24694339e+02,  8.62095906e+01],
        [ 1.69471687e-01,  1.99983894e+01,  8.62095906e+01,  2.16130608e+03],
    ])
    n_run4 = 2295198

    cov_run5 = np.array([
        [ 3.18657505e-01, -7.93943144e-03,  2.70711413e+00,  1.27049737e-01],
        [-7.93943144e-03,  5.81624810e+00, -5.50813585e+00,  1.95119974e+01],
        [ 2.70711413e+00, -5.50813585e+00,  9.24057080e+02,  8.29724938e+01],
        [ 1.27049737e-01,  1.95119974e+01,  8.29724938e+01,  2.14144404e+03]
    ])
    n_run5 = 2634946

    cov_list = [cov_run1, cov_run2, cov_run3, cov_run4, cov_run5]
    n_list = np.array([n_run1, n_run2, n_run3, n_run4, n_run5])

    weights = n_list / n_list.sum()

    cov_avg = np.zeros_like(cov_list[0])
    for w, cov in zip(weights, cov_list):
        cov_avg += w * cov

    # conversione a unità SI: x,y in mm->m, px,py in µrad->rad
    scale = np.array([1e-3, 1e-3, 1e-6, 1e-6])
    cov_avg = cov_avg * np.outer(scale, scale)

    eigvals = np.linalg.eigvalsh(cov_avg)
    # print("min eigenvalue:", eigvals.min())
    # print("\nCovariance matrix (averaged):")
    # print(cov_avg)
    return cov_avg, eigvals, n_list


def mean_vector():
    mu_run1 = np.array([-0.12e-3, 0.64e-3, -0.93e-6, 6.79e-6])  # [mu_x, mu_y, mu_px, mu_py] del run 1
    mu_run2 = np.array([-0.25e-3, 0.96e-3, -1.61e-6, 8.18e-6])
    mu_run3 = np.array([-0.40e-3, 0.58e-3, -2.79e-6, 6.05e-6])
    mu_run4 = np.array([-0.32e-3, -0.01e-3, -1.74e-6, -1.03e-6])
    mu_run5 = np.array([-0.29e-3, -0.03e-3, -1.38e-6, -0.03e-6])

    mu_list = [mu_run1, mu_run2, mu_run3, mu_run4, mu_run5]
    n_list = np.array([4263962, 3807353, 2651593, 2295198, 2634946])
    weights = n_list / n_list.sum()

    mu_avg = np.zeros_like(mu_list[0])
    for w, mu in zip(weights, mu_list):
        mu_avg += w * mu

    return mu_avg