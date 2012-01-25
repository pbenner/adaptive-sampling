function plot_moments(results, varargin)
% plot_moments(results [, parameters])
%
% plot moments of the MPM posterior
%
% results: structure of results returned by
% algorithm specified in
%   Poppe, Benner & Elze (2012)
%
% the following parameters are implemented:
%
% parameter, default [: description]
%  'marginal', 1: plot full marginal distribution
%  'marginal_colormap', flipud(gray()): colormap for plotting marginals
%  'marginal_colorbar_location', 'NorthOutside': colorbar for marginals;
%    see help colorbar for more details about locations;
%    'marginal_colorbar_location', 'off' switches colorbar off
%  'autoclip_marginal', 1: clip y axis of the marginals to sensible values
%    if not set, the whole range between 0 and 1 is plotted
%  'moment1_plotprops', 'r': plot properties for first moment
%  'moment2_plotprops', 'r--': plot properties for second moment
%  'std', 1: plot standard deviation (sqrt of 2nd moment)
%    around the first moment (expected value)
%  'bprob', 0: plot break probabilities
%  'bprob_color', 'g': line color for break probabilities
%
% Examples:
%  counts = [10 11 10 10 12 10 11 20 21 19 19 20 20 19 21 20; ...
%    90 89 90 90 88 90 89 80 79 81 81 80 80 81 79 80];
%  r = adsamp(counts, 'marginal_step', 0.005)
%  plot_moments(r);
%  plot_moments(r, 'bprob', 1, 'marginal_colorbar_location', 'off');

% options:
p = inputParser;
p.addParamValue('marginal', 1, @isscalar);
p.addParamValue('marginal_colormap', flipud(gray()), @isnumeric);
p.addParamValue('marginal_colorbar_location', 'NorthOutside', @ischar);
p.addParamValue('autoclip_marginal', 1, @isscalar);
p.addParamValue('moment1_plotprops', 'r', @ischar);
p.addParamValue('moment2_plotprops', 'r--', @ischar);
p.addParamValue('std', 1, @isscalar);
p.addParamValue('bprob', 0, @isscalar);
p.addParamValue('bprob_color', 'g', @ischar);
p.KeepUnmatched = true;
p.parse(varargin{:});

if ~isfield(results, 'moments')
	error('results structure does not contain moments')
end
m1 = results.moments(1, :);
if length(results.moments(:,1)) > 1
	has_moment2 = 1;
	variance = results.moments(2, :) - results.moments(1, :).^2;
	stdev = sqrt(variance);
	std_top = m1 + stdev;
	std_bot = m1 - stdev;
end

if p.Results.marginal
	if ~isfield(results, 'marginals')
		error('results structure does not contain marginals')
	end
	imagesc(...
		[1 length(results.marginals(:,1))], ...
		[1 0], ...
		flipud(results.marginals'));
	set(gca, 'YDir', 'normal');
	colormap(p.Results.marginal_colormap);
	if ~strcmp(p.Results.marginal_colorbar_location, 'off')
		colorbar('location', p.Results.marginal_colorbar_location)
	end
	xlim([1 length(results.marginals(:,1))]);
	if p.Results.autoclip_marginal
		if ~has_moment2
			warning('cannot autoclip marginals: results structure does not contain 2nd moment')
		end
		ylimtop = min(max(m1 + 1.5*stdev), 1);
		ylimbot = max(min(m1 - 1.5*stdev), 0);
		ylim([ylimbot ylimtop]);
	end
	
	hold on
end

% plot expected value (1st moment):
plot(m1, p.Results.moment1_plotprops);

if p.Results.std
	if ~has_moment2
		error('results structure does not contain 2nd moment')
	end
	hold on
	plot(std_top, p.Results.moment2_plotprops);
	plot(std_bot, p.Results.moment2_plotprops);
end

if p.Results.bprob
	if ~isfield(results, 'bprob')
		error('results struct does not contain bprob')
	end
	hold on
	ax1 = gca;
	col = p.Results.bprob_color;
	ax2 = axes('Position',get(ax1,'Position'),...
		'XAxisLocation','top',...
		'xticklabel',[], ...
		'YAxisLocation','right',...
		'Color','none',...
		'XColor', col, 'YColor', col);
	bprob = results.bprob;
	bprob(1) = nan;
	bprob(end) = nan;
	line(1:length(bprob), bprob, 'Color', col, 'Parent', ax2);
	xlim([1 length(bprob)]);
	ylimits = get(ax1,'YLim');
	yinc = (ylimits(2)-ylimits(1))/5;
	set(ax1,'YTick',[ylimits(1):yinc:ylimits(2)]);
	ylimits = get(ax2,'YLim');
	yinc = (ylimits(2)-ylimits(1))/5;
	set(ax2,'YTick',[ylimits(1):yinc:ylimits(2)]);
end

hold off

