function result = default_beta(L)
%DEFAULT_BETA - Computes a set of default values for the beta parameter on log scale

%
  result = ones(1,L)/L;
  for i = 1:L
      result(i) = log(result(i)) - lnchoose(L-1,i-1);
  end
%
end % default_beta
