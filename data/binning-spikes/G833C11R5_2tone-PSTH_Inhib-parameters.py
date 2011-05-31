
def parameters(K, L):
    alpha_v = np.ones([K, L])
    alpha = generate_alpha(alpha_v)
    beta  = generate_beta(range(1,26), L)
    gamma = generate_gamma(L)

    return alpha, beta, gamma
