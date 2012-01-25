function plot_binning(results, varargin)
% plot_binning(results [, parameters])
%
% plot results returned by algorithm specified in
%   Poppe, Benner & Elze (2012)
%
% The following parameters are implemented:
% all parameters of function plot_moments, and in addition:
%
% parameter, default [: description]
%  'bar_color', [0 0.5 0]: color of the bars of the class posteriors
%
% Examples:
%  counts = [1 2 2 3 2 4 4 5 3 4 6 5 4 5 4 5 5 4 3 4 4 5 5 6 6 6 7 8 8 8 7 5 4 4 2 2; ...
%    9 8 8 7 8 6 6 5 7 6 4 5 6 5 6 5 5 6 7 6 6 5 5 4 4 4 3 2 2 2 3 5 6 6 8 8];
%  r = adsamp(counts)
%  plot_binning(r);
%  plot_binning(r, 'bprob', 1);

p = inputParser;
p.addParamValue('bar_color', [0 0.5 0], @(x) ischar(x) || length(x) == 3);
p.KeepUnmatched = true;
p.parse(varargin{:});

figure('Position', [10 10 580 800]);

axes('Position', [0.07 0.51 0.86 0.45]);
plot_moments(results, 'XAxisLocation', 'top', varargin{:});

axes('Position', [0.07 0.05 0.86 0.44]);
if ~isfield(results, 'mpost')
	error('results structure does not contain mpost')
end
bar(results.mpost, 'FaceColor', p.Results.bar_color);
xlim([1 length(results.mpost)]);
grid on
