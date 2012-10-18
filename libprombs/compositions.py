def ruleGen(n, m, sigma):
    """
    Generates all interpart restricted compositions of n with first part
    >= m using restriction function sigma. See Kelleher 2006, 'Encoding
    partitions as ascending compositions' chapters 3 and 4 for details.
    """
    a = [0 for i in range(n + 1)]
    k = 1
    a[0] = m - 1
    a[1] = n - m + 1
    while k != 0:
        x = a[k - 1] + 1
        y = a[k] - 1
        k -= 1
        while sigma(x) <= y:
            a[k] = x
            x = sigma(x)
            y -= x
            k += 1
        a[k] = x + y
        yield a[:k + 1]

ruleGen(n, 1, lambda x: 1)
