function result = prombsTest(m)
%
  if nargin < 1
    m = 5;
  end

  g = [1 2 3 4 5]
  f = [[1 2 3 4 5]; [0 1 2 3 4]; [0 0 1 2 3]; [0 0 0 1 2]; [0 0 0 0 1]]

  % correct answer is [5 40 63 32 5]
  %
  result = exp(prombs(log(g), triu(log(f)), m));

end % prombsTest
