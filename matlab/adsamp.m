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
%  'kl_psi', 1
%  'kl_multibin', 0
%  'effective_counts', 0
%  'effective_posterior_counts', 0
%  'density', 1: compute full density distribution
%  'density_step', 0.01: step size for the density distribution
%  'density_range', [0.0 1.0]: limit range for the density distribution
%  'epsilon', 0.00001: precision for the extended prombs
%  'threads', getNumberOfCores: number of threads
%  'stacksize', 256*1024: thread stack size
%  'algorithm', 'prombs': select an algorithm 
%      [prombs | mgs]
%  'which', 0: for which event to compute the binning
%  'hmm', 0: use hidden Markov model
%  'rho', 0.4: cohesion parameter for the hidden Markov model
%  'samples', [100 2000]
%
%
% Example:
%  counts = [ 1 1 2 1; 1 1 1 1 ];
%  adsamp(counts)
%  r = adsamp(counts, 'density_range', [0.0 0.5]);
%  r.density

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
p.addParamValue('kl_psi', 1, @isscalar);
p.addParamValue('kl_multibin', 0, @isscalar);
p.addParamValue('effective_counts', 0, @isscalar);
p.addParamValue('effective_posterior_counts', 0, @isscalar);
p.addParamValue('density', 1, @isscalar);
p.addParamValue('density_step', 0.01, ispos);
p.addParamValue('density_range', [0.0 1.0], isinterval);
p.addParamValue('epsilon', 0.00001, ispos);
p.addParamValue('threads', getNumberOfCores, @isscalar);
p.addParamValue('stacksize', 256*1024, ispos);
p.addParamValue('algorithm', 'prombs', @ischar);
p.addParamValue('which', 0, @isscalar);
p.addParamValue('hmm', 0, @isscalar);
p.addParamValue('rho', 0.4, @isscalar);
p.addParamValue('samples', [100 2000], ispair);
p.KeepUnmatched = true;
p.parse(varargin{:});

options = p.Results;
switch p.Results.algorithm
	case 'prombs'
		options.algorithm = 0;
	case 'mgs'
		if p.Results.utility
			error('utility cannot be calculated for algorithm mgs')
		end
		options.algorithm = 1;
	otherwise
		error(['algorithm ' p.Results.algorithm ' is not implemented'])
end

result  = samplingUtility(countstat, args.Results.alpha, args.Results.beta, args.Results.gamma, options);
