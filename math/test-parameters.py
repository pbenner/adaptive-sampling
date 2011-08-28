
def parameters(K, L):
    alpha_v = [ L * [ 1 ], L * [ 2 ] ]
    alpha = generate_alpha(alpha_v)
    beta  = generate_beta([], L)
    gamma = generate_gamma(L)

    return alpha, beta, gamma
