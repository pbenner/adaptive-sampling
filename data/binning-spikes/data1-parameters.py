
def parameters(K, L):
    alpha_v = np.ones([K, L])
    alpha_v[1] *= 32
    alpha = generate_alpha(alpha_v)
    beta  = generate_beta(range(1,20+1), L)
    gamma = generate_gamma(L)

    return alpha, beta, gamma
