import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
style.use('ggplot')


def finite_difference_solver(M, N):
    sigma = 0.40
    T = 1.00
    r = 0.05

    a = r - (sigma**2)/2
    b = (sigma**2)/2
    c = -r

    x_min = np.log(50)
    x_max = np.log(100)
    s_min = 0
    s_max = T
    x_range, dx = np.linspace(x_min, x_max, N, retstep=True)
    s_range, ds = np.linspace(s_min, s_max, M, retstep=True)

    V = np.zeros((M,N))
    V[0,:] = (np.exp(x_range) - 50) * (100 - np.exp(x_range))
    for m_ix, m in enumerate(s_range):
        if m_ix == len(s_range) - 1:
            continue
        for n_ix, n in enumerate(x_range):
            if n_ix == 0:
                V[m_ix, n_ix] = 0 # lower boundary condition
                continue
            if n_ix == len(x_range) - 1:
                V[m_ix, n_ix] = 0 # upper boundary condition
                continue

            prev_left = V[m_ix, n_ix - 1]
            prev_right = V[m_ix , n_ix + 1]
            prev_center = V[m_ix, n_ix]

            V[m_ix+1, n_ix] = ds*a*(prev_right - prev_left)/(2*dx) + ds*b*(prev_right - 2*prev_center + prev_left)/(dx**2) + prev_center*(ds*c + 1)

    return V, x_range, s_range



if __name__ == '__main__':
    N = 100
    k = 2
    x_range, dx = np.linspace(np.log(50), np.log(100), N, retstep=True)
    ds = k * dx**2
    M = int(round(1/ds))
    V, x_range, s_range = finite_difference_solver(M, N)
