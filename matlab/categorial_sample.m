function x = categorial_sample(p)
%
% returns a sample from the categorial distribution with
% parameters p
%
%   - p: vector of length n with probabilities for each
%        event i = 1, 2, ..., n
%

cdf = cumsum(p(:));
x   = sum(cdf < rand*cdf(end)) + 1;
