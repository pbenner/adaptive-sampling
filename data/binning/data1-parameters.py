
def parameters(K, L):
    alpha = generate_alpha(np.ones([K, L]))
    beta  = generate_beta([], L)
    gamma = generate_gamma(L)

    return alpha, beta, gamma
