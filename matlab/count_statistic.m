function out = count_statistic(in)
% DEFAULT_OPTIONS
%   Compute the full count statistic.
%   TODO: This is a naive implementation! Use the outer product!
%
  [K, L] = size(in);
  out    = [];

  for k = 1:K
    m = compute_counts(in(k,:));
    out = [out; m];
  end

end % count_statistic

function out = compute_counts(in)
  [K, L] = size(in);

  out = zeros(L, L, 'double');

  for i = 1:L
    for j = i:L
      sum = 0.0;
      for k = i:j
        sum = sum + in(k);
      end
      out(i,j) = sum;
    end
  end

end % compute_counts
