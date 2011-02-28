def preplot(result):
    samples = len(result['samples'])
    result['samples'] = np.array(result['samples']) - 17
    x = np.arange(-17, 18, 1)
    title = "(b) "+str(samples)+" Samples"
    return x, title

def postplot(ax, p, result, options):
    [ax11, ax12, ax21, ax22, ax31, ax32] = ax 
    [ p11,  p12,  p21,  p22,  p31,  p32] = p
    ax11.set_xlim(-17,17)
    ax11.set_xlabel('deviant offset', font)
    ax11.set_ylabel(r'$P(s_i|{\bf X}_n, {\bf Y}_n)$', font)
    ax11.set_ylim(0,1)
    ax11.legend([p11, p12], [r'$P(s_i|{\bf X}_n, {\bf Y}_n)$', 'Ground Truth'], loc=4, prop=smallfont)
    ax12.set_ylabel('Ground Truth', font)
    ax21.set_xlim(-17,17)
    ax21.set_ylabel('Counts', font)
    ax21.set_xlabel('deviant offset',  font)
    ax21.legend([p22], ['$U_x$'], loc=4, prop=smallfont)
    ax22.set_ylabel(r'$U_x$', font)
    ax31.set_ylabel(r'$P(m_B|{\bf X}_n, {\bf Y}_n)$', font)
    ax31.set_xlabel(r'$m_B$',  font)

def utilitypreplot(result):
    x = np.array(range(len(result['samples'])-len(result['utility']), len(result['samples'])))
    y = np.array(range(0, len(result['utility'][0]))) - 17
    return x, y, ''

def utilitypostplot(ax, p, cb):
    ax.set_xlabel('i')
    ax.set_ylabel('deviant offset')
    cb.set_label(r'normalized $U_x$')
