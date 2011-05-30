
def parameters(K, L):
    alpha_v = np.ones([K, L])
    alpha_v[0][-1] = 100
    alpha_v[1][ 0] = 100
    alpha = generate_alpha(alpha_v)
    beta  = generate_beta([], L)
    gamma = generate_gamma(L)

    return alpha, beta, gamma
