def preplot(result, options):
    x = np.arange(0, len(result['moments'][0]), 1)
    return x, ''

def postplot(ax, p, result, options):
    [ax11, ax12, ax21, ax22] = ax 
    [ p11,  p12,  p21,  p22] = p
    if result['bprob'] and options['bprob']:
        ax12.set_ylabel(r'$P(\Rsh_i|E)$', font)
    ax11.set_ylabel(r'$P(s_i|E)$', font)
    ax11.set_xlabel(r'$x$',  font)
    ax21.set_ylabel(r'$P(M = m_B|E)$', font)
    ax21.set_xlabel(r'$m_B$',  font)
