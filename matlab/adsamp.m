function result = adsamp(counts, varargin)
% result = adsamp(counts [, parameters])
%
% algorithm specified in
%   Poppe, Benner & Elze (2012)
%
% counts: matrix of counts; each line one response option
%  dimensions: K rows, L columns
%
% the following parameters are implemented:
%
% parameter, default [: description]
%  'alpha', default_alpha(ones(K, L)): "pseudo counts"
%  'beta', ones(1, L): relative class weights
%  'gamma', default_gamma(L): a priori importance of each consecutive bin
%  'n_moments', 2: compute the first N>=2 moments
%  'model_posterior', 2
%  'bprob', 1
%  'utility', 1
%  'differential_entropy', 1
%  'multibin_entropy', 0
%  'effective_counts', 0
%  'marginal', 1: compute full marginal distribution
%  'marginal_step', 0.01: step size for the marginal distribution
%  'marginal_range', [0.0 1.0]: limit range for the marginal distribution
%  'epsilon', 0.00001: epsilon for entropy estimations
%  'threads', 1: number of threads
%  'stacksize', 256*1024: thread stack size
%  'algorithm', 'prombs': select an algorithm 
%      [mgs | prombstree | prombs]
%  'which', 0: for which event to compute the binning
%  'samples', [100 2000]
%
%
% Example:
%  counts = [ 1 1 2 1; 1 1 1 1 ];
%  adsamp(counts)
%  r = adsamp(counts, 'marginal_range', [0.0 0.5]);
%  r.marginals

[K, L] = size(counts);
countstat  = count_statistic(counts);

ispair = @(x) isnumeric(x) && length(x) == 2;
isinterval = @(x) ispair(x) && 0 <= x(1) && x(1)<x(2) && x2<=1;
ispos = @(x) isscalar(x) && x>0;

args = inputParser;
args.addParamValue(...
	'alpha', ...
	default_alpha(ones(K, L)), ...
	@(x) isnumeric(x) && all(size(x) == [K*L L]));
args.addParamValue('beta', ones(1, L), @(x) isvector(x) && length(x)==L);
args.addParamValue(...
	'gamma', ...
	default_gamma(L), ...
	@(x) isnumeric(x) && all(size(x) == L));
args.KeepUnmatched = true;
args.parse(varargin{:});

% options:
p = inputParser;
p.addParamValue('n_moments', 2, @isscalar);
p.addParamValue('model_posterior', 2, @isscalar);
p.addParamValue('bprob', 1, @isscalar);
p.addParamValue('utility', 1, @isscalar);
p.addParamValue('differential_entropy', 1, @isscalar);
p.addParamValue('multibin_entropy', 0, @isscalar);
p.addParamValue('effective_counts', 0, @isscalar);
p.addParamValue('marginal', 1, @isscalar);
p.addParamValue('marginal_step', 0.01, ispos);
p.addParamValue('marginal_range', [0.0 1.0], isinterval);
p.addParamValue('epsilon', 0.00001, ispos);
p.addParamValue('threads', 1, @isscalar);
p.addParamValue('stacksize', 256*1024, ispos);
p.addParamValue('algorithm', 'prombs', @ischar);
p.addParamValue('which', 0, @isscalar);
p.addParamValue('samples', [100 2000], ispair);
p.KeepUnmatched = true;
p.parse(varargin{:});

options = p.Results;
switch p.Results.algorithm
	case 'prombs'
		options.algorithm = 0;
	case 'prombstree'
		options.algorithm = 1;
	case 'mgs'
		if p.Results.utility
			error('utility cannot be calculated for algorithm msg')
		end
		options.algorithm = 2;
	otherwise
		error(['algorithm ' p.Results.algorithm ' is not implemented'])
end


result  = adaptive_sampling(countstat, args.Results.alpha, args.Results.beta, args.Results.gamma, options);

