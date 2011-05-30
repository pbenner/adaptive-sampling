def preplot(x, timings, result, options):
    return x, 1, timings, ''

def postplot(ax, p, result, options):
    [ax11, ax12, ax21, ax22, ax31, ax32] = ax 
    [ p11,  p12,  p21,  p22,  p31,  p32] = p
    ax11.set_xlabel(r'$t$')
    ax11.set_ylabel('Trial')
    if result['bprob'] and options['bprob']:
        ax12.set_ylabel(r'$P(\Rsh_i|E)$', font)
    ax21.set_ylabel(r'$P(s_i|E)$', font)
    ax21.set_xlabel(r'$t$',  font)
    ax31.set_ylabel(r'$P(m_B|E)$', font)
    ax31.set_xlabel(r'$m_B$',  font)
    ax31.set_xlim(0,30)
