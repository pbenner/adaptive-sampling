
def parameters(K, L):
    alpha_v = np.ones([K, L])
    alpha_v[0] *= 32
    alpha = generate_alpha(alpha_v)
    beta  = generate_beta([3,4], L)
    gamma = generate_gamma(L)

    return alpha, beta, gamma
