def preplot(result, options):
    x = np.arange(1, 36, 1)
    samples = len(result['samples'])
    result['samples'] = np.array(result['samples']) + 1
    if samples == 0:
        title = "(a) "+str(samples)+" Samples"
    elif samples == 10:
        title = "(b) "+str(samples)+" Samples"
    elif samples == 100:
        title = "(c) "+str(samples)+" Samples"
    elif samples == 200:
        title = "(d) "+str(samples)+" Samples"
    elif samples == 300:
        title = "(e) "+str(samples)+" Samples"
    else:
        title = "(f) "+str(samples)+" Samples"
    return x, title

def postplot(ax, p, result, options):
    import re
    font = {'family'     : 'serif',
            'weight'     : 'normal',
            'size'       : 16 }
    [ax11, ax12, ax21, ax22, ax31, ax32] = ax
    [ p11,  p12,  p21,  p22,  p31,  p32] = p
    # adjust style
    ax11.lines[2].set_linewidth(2.5)
    ax12.lines[0].set_linewidth(0.0)
    ax12.lines[0].set_linestyle('--')
    ax12.lines[0].set_marker('_')
    ax12.lines[0].set_markersize(12)
    ax12.lines[0].set_mew(2.5)
    # set limits and add legends
    ax11.set_xlim(1,35)
    ax11.set_xlabel(r'$x$', font)
    ax11.set_ylabel(r'$p_{x,s}$', font)
    ax11.set_ylim(0,1)
    ax11.legend([p11, p12], [r'$\mathbb{E}( p_{x,s} \mid {\bf Y}^\mathcal{X}_{{\bf n}_\mathcal{X}})$', 'Ground truth'], ncol=2, mode='expand', loc=3, frameon=False, borderaxespad=0., prop=font)
    ax12.set_ylabel('Ground truth', font)
    ax21.set_xlim(1,35)
    ax21.set_ylabel('Counts', font)
    ax21.set_xlabel(r'$x$',  font)
    ax21.legend([p22], ['$U(x)$'], loc=4, prop=font)
    ax22.set_ylabel(r'$U(x)$', font)
    if ax31:
        ax31.set_ylabel(r'$P(M = m \mid {\bf Y}^\mathcal{X}_{{\bf n}_\mathcal{X}})$', font)
        ax31.set_xlabel(r'$m$',  font)
        ax32.set_axis_off()
