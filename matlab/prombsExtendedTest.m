function result = prombsExtendedTest(m)
%
  if nargin < 1
    m = 5;
  end

  g = [1 2 3 4 5]
  f = [[1 2 3 4 5]; [0 1 2 3 4]; [0 0 1 2 3]; [0 0 0 1 2]; [0 0 0 0 1]]
  h = [[1 2 3 4 5]; [0 1 2 3 4]; [0 0 1 2 3]; [0 0 0 1 2]; [0 0 0 0 1]]

  % correct answer is [25.0063 200.0500 315.0788 160.0400 25.0063]
  %
  result = exp(prombsExtended(log(g), triu(log(f)), h, 0.0001, m));

end % prombsTest
