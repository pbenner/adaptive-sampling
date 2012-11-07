function result = example3()
%
  bare_counts = [29 6 5 6 7 9 16 22 29 20 12 9 10 18 9 4 3 3 3 6 10 16 19 21 24 17 22 11 11 7 6 6 5 6 27;
	              0 0 0 0 0 0 1 2 6 8 5 6 6 11 20 8 6 7 7 13 14 13 12 8 6 4 2 0 1 0 0 0 0 0 0];

  % K: responses
  % L: stimuli
  [K, L]  = size(bare_counts);

  counts  = count_statistic(bare_counts);
  alpha   = default_alpha(ones(K, L));
  beta    = default_beta(L);
  gamma   = default_gamma(L);

  options = default_options();

  result  = samplingUtility(counts, alpha, beta, gamma, options);

end % example3
