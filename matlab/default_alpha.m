function out = default_alpha(in)
% DEFAULT_OPTIONS
%   Convert a matrix to a list of matrices with default alpha values.
%
  [K, L] = size(in);
  out    = [];

  for k = 1:K
    m = generate_alpha(in(k,:));
    out = [out; m];
  end

end % default_alpha

function out = generate_alpha(in)
  [K, L] = size(in);

  out = zeros(L, L, 'double');

  for i = 1:L
    for j = i:L
      sum = 0.0;
      for k = i:j
        sum = sum + in(k);
      end
      out(i,j) = sum/double(j-i+1);
    end
  end

end % generate_alpha
