function result = demo()
%
  bare_counts =[0,3,5,8,9,5,7,9,9,9,9,11,12,4,3,0,0,2,8,14,14,8,7,6,6,6,5,7,5,3,0;
                0,0,0,1,2,0,0,1,2,3,3,4,7,6,2,0,1,4,8,5,3,1,1,0,1,0,0,1,0,0,0;]

  % K: responses
  % L: stimuli
  [K, L]  = [2, 31]

  counts  = count_statistic(doublebare_counts);
  alpha_v = ones(K, L);
  alpha_v(1,1) = 100;
  alpha_v(1,L) = 100;
  alpha   = default_alpha(alpha_v);
  % Alpha enthaelt jetzt:
  % squeeze(alpha(1,:,:)) Gradient zum Rand
  % squeeze(alpha(2,:,:)) Uniform
  % damit ist der Prior fuer den Rand gesetzt.
  % Jetzt noch die Mitte, jedes bin welches die Mitte enthaelt
  % soll auf 1/3 gesetzt werden:
  alpha(1,1:16,16:end) = 100;
  alpha(2,1:16,16:end) = 200;
  beta    = ones(1, L);
  gamma   = default_gamma(L);

  options = default_options();

  result  = samplingUtility(counts, alpha, beta, gamma, options);

end % demo
