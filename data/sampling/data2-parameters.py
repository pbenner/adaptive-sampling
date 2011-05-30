
def parameters(K, L):
    alpha_v = np.ones([K, L])
    alpha = generate_alpha(alpha_v)
    beta  = generate_beta([1,2,3,4,5], L)
    gamma = generate_gamma(L)

    return alpha, beta, gamma
