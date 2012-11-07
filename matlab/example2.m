function result = example2()
%
  bare_counts = [6 1 2 1 1 2 1 2 2 1 2 2 6 2 0 0 1 1 1 0 3 0 2 5 3 3 2 2 1 2 1 1 2 2 6;
	        	 0 0 0 0 0 0 0 0 1 0 0 0 1 5 3 2 2 1 1 2 2 3 5 1 1 1 0 0 0 0 0 0 0 0 0];

  % K: responses
  % L: stimuli
  [K, L]  = size(bare_counts);

  counts  = count_statistic(bare_counts);
  alpha   = default_alpha(ones(K, L));
  beta    = default_beta(L);
  gamma   = default_gamma(L);

  options = default_options();

  result  = samplingUtility(counts, alpha, beta, gamma, options);

end % example2
