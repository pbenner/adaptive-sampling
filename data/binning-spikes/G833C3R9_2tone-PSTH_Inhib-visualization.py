def preplot(x, timings, result, options):
    x = map(lambda y: y/100.0, x)
    timings = map(lambda x: map(lambda y: y/100.0, x), timings)
    return x, 0.1*0.001, timings, ''

def postplot(ax, p, result, options):
    [ax11, ax12, ax21, ax22, ax31, ax32] = ax 
    [ p11,  p12,  p21,  p22,  p31,  p32] = p
    ax11.lines[0].set_marker(',')
    ax21.lines[0].set_linestyle('-')
    ax21.lines[1].set_linestyle('-')
    ax21.lines[0].set_linewidth(0.5)
    ax21.lines[1].set_linewidth(0.5)
    ax21.lines[2].set_linewidth(1.5)
    ax11.set_xlim(0,300)
    ax21.set_xlim(0,300)
    ax11.set_xlabel(r'$t$ [ms]')
    ax11.set_ylabel('Trial')
    if result['bprob'] and options['bprob']:
        ax12.set_ylabel(r'$P(\Rsh_t|E)$', font)
    else:
        ax12.set_axis_off()
    ax21.set_ylabel(r'$r(t) = P(s_t|E)/dt$ [pps]', font)
    ax21.set_xlabel(r'$t$ [ms]',  font)
    ax31.set_ylabel(r'$P(m_B|E)$', font)
    ax31.set_xlabel(r'$m_B$',  font)
    ax31.set_xlim(0,30)
    ax22.set_axis_off()
    ax32.set_axis_off()
