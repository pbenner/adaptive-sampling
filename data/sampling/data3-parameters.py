
def parameters(K, L):
    alpha_v = np.ones([K, L])
    alpha = generate_alpha(alpha_v)
    beta  = generate_beta([], L)
    gamma = generate_gamma(L)

    return alpha, beta, gamma
