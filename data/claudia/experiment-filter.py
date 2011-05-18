def gainFilter(gain, counts, result):
    for i in range(0, len(gain)):
        if counts[i] > 15:
            print "stopping to sample at position "+str(i)
            gain[i] = 0
    return gain
