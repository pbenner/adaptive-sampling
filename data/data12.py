def preplot(result):
    x = np.arange(0, 31, 1)
    title = str(len(result['samples']))+" Samples"
    return x, title

def postplot(ax, p, result, options):
    [ax11, ax12, ax21, ax22, ax31, ax32] = ax 
    [ p11,  p12,  p21,  p22,  p31,  p32] = p
    ax11.set_xlabel('t', font)
    ax11.set_ylabel(r'$P(s_i|E)$', font)
    ax11.set_ylim(0,1)
    ax11.legend([p11, p12], ['$P(s_i|E)$', 'Ground Truth'], loc=4, prop=smallfont)
    ax12.set_ylabel('Ground Truth', font)
    ax21.set_ylabel('Counts', font)
    ax21.set_xlabel('input',  font)
    ax21.legend([p22], ['$U_x$'], loc=4, prop=smallfont)
    ax22.set_ylabel(r'$U_x$', font)
    ax31.set_ylabel(r'$P(m_B|E)$', font)
    ax31.set_xlabel(r'$m_B$',  font)
    ax31.set_xlim(0,30)


def utilitypreplot(result):
    x = range(1, len(result['utility'])+1)
    y = np.array(range(0, len(result['utility'][0])))
    return x, y, ''

def utilitypostplot(ax, p, cb):
    ax.set_xlabel('i')
    ax.set_ylabel('input')
    cb.set_label(r'$U_x$')
