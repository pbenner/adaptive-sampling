function result = example4()
%
  bare_counts = [0,3,5,8,9,5,7,9,9,9,9,11,12,4,3,0,0,2,8,14,14,8,7,6,6,6,5,7,5,3,0;
                 0,0,0,1,2,0,0,1,2,3,3,4,7,6,2,0,1,4,8,5,3,1,1,0,1,0,0,1,0,0,0];

  % K: responses
  % L: stimuli
  [K, L]  = size(bare_counts);

  counts  = count_statistic(bare_counts);
  alpha_v = ones(K, L);
  alpha_v(1,1) = 100;
  alpha_v(1,L) = 100;
  alpha   = default_alpha(alpha_v);
  % Alpha contains:
  % squeeze(alpha(1,:,:)): gradient to the sides
  % squeeze(alpha(2,:,:)): uniform
  % The center should be 1/3 chance:
  alpha(1,1:16,16:end) = 100;
  alpha(2,1:16,16:end) = 200;
  beta    = default_beta(L);
  gamma   = default_gamma(L);

  options = default_options();

  result  = samplingUtility(counts, alpha, beta, gamma, options);

end % example4
