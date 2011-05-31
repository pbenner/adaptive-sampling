
def parameters(K, L):
    alpha_v = np.ones([K, L])
    alpha_v[1] *= 32
    alpha = generate_alpha(alpha_v)
    beta  = generate_beta([2,3], L)
    gamma = generate_gamma(L)

    return alpha, beta, gamma
