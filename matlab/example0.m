function result = example0()
%
  bare_counts = [ 1 1 2 1;
                  1 1 1 1 ];

  % K: responses
  % L: stimuli
  [K, L]  = size(bare_counts);

  counts  = count_statistic(bare_counts);
  alpha   = default_alpha(ones(K, L));
  beta    = default_beta(L);
  gamma   = default_gamma(L);

  options = default_options();

  result  = samplingUtility(counts, alpha, beta, gamma, options);

end % example0
