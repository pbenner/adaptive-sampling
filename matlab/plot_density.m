function plot_density(results, varargin)
% plot_density(results [, parameters])
%
% plot MPM density
%
% results: structure of results returned by algorithm specified in
%   Poppe, Benner & Elze (2012)
%
% The following parameters are implemented:
%
% parameter, default [: description]
%  'density', 1: plot full density distribution (if available);
%    if no density available, the first two moments are plotted
%  'density_colormap', flipud(gray(256)): colormap for plotting densitys
%  'density_colorbar_location', 'off': colorbar for densitys;
%    see help colorbar for more details about locations;
%    'density_colorbar_location', 'off' switches colorbar off
%  'autoclip_density', 1: clip y axis of the densitys to sensible values
%    if not set, the whole range between 0 and 1 is plotted
%  'normalize_density_plot', 0: normalize densitys so that different
%    plots are easier to compare; if z are densitys, normalized densitys
%    n are n = log(1+z) && n <= 5
%  'median_plotprobs', 'r': plot properties of median (if density 
%    available)
%  'quartiles', 1: plot 1st and 3rd quartile
%  'quartiles_plotprobs', 'r--': plot properties of quartiles (apart from
%    median)
%  'outer_quantiles', 1: plot quantiles 0.025 and 0.975
%  'outer_quantiles_plotprobs', 'r:': plot properties of quantiles
%    0.025 and 0.975
%  'add_moments', 0: add first two moments to the plot 
%  'moment1_plotprops', 'r': plot properties of first moment
%  'moment2_plotprops', 'r--': plot properties of second moment
%  'std', 1: plot standard deviation (sqrt of 2nd moment)
%    around the first moment (expected value), if moment plotting
%    is active
%  'bprob', 0: plot break probabilities
%  'bprob_color', [0 0.5 0]: line color for break probabilities
%  'show_legend', 1: show legend
%  'legend_location', 'Best': location of the legend (see help legend)
%  'XAxisLocation', 'bottom'
%  'X' : X Labels
%
% Examples:
%  counts = [10 11 10 10 12 10 11 20 21 19 19 20 20 19 21 20; ...
%    90 89 90 90 88 90 89 80 79 81 81 80 80 81 79 80];
%  r = adsamp(counts, 'density_step', 0.001)
%  plot_density(r);
%  plot_density(r, ...
%    'density_colorbar_location', 'NorthOutside', ...
%    'density_colormap', winter);
%
%  counts = [1 2 2 3 2 4 4 5 3 4 6 5 4 5 4 5 5 4 3 4 4 5 5 6 6 6 7 8 8 8 7 5 4 4 2 2; ...
%    9 8 8 7 8 6 6 5 7 6 4 5 6 5 6 5 5 6 7 6 6 5 5 4 4 4 3 2 2 2 3 5 6 6 8 8];
%  r = adsamp(counts)
%  plot_density(r, 'bprob', 1);

% options:
p = inputParser;
p.addParamValue('density', 1, @isscalar);
p.addParamValue('density_colormap', flipud(gray(256)), @isnumeric);
p.addParamValue('density_colorbar_location', 'off', @ischar);
p.addParamValue('autoclip_density', 1, @isscalar);
p.addParamValue('normalize_density_plot', 0, @isscalar);
p.addParamValue('median_plotprobs', 'r', @ischar);
p.addParamValue('quartiles', 1, @isscalar);
p.addParamValue('quartiles_plotprobs', 'r--', @ischar);
p.addParamValue('outer_quantiles', 1, @isscalar);
p.addParamValue('outer_quantiles_plotprobs', 'r:', @ischar);
p.addParamValue('add_moments', 0, @isscalar);
p.addParamValue('moment1_plotprops', 'r', @ischar);
p.addParamValue('moment2_plotprops', 'r--', @ischar);
p.addParamValue('std', 1, @isscalar);
p.addParamValue('bprob', 0, @isscalar);
p.addParamValue('bprob_color', [0 0.5 0], @(x) ischar(x) || length(x) == 3);
p.addParamValue('show_legend', 1, @isscalar);
p.addParamValue('legend_location', 'Best', @(x) ischar(x) || length(x) == 4);
p.addParamValue('XAxisLocation', 'bottom', @ischar);
p.addParamValue('X',[1:length(results.bprob)],@(x) isnumeric(x) && length(x) == length(results.bprob));
p.addParamValue('YLim',[0,1],@isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});

% initialize legend handles:
legend_handles = [];
legend_strings = {};

if isfield(results, 'density')
	% calculate quantiles [.025 .25 .50 .75 .975]:
	n = length(results.density(:, 1));
	quantiles_p = [.025 .25 .50 .75 .975];
	quantiles = zeros(5, n);
	for i=1:n
		v = results.density(i, :);
		% normalized cumsum:
		ncumsum = cumsum(v)/sum(v);
		% make data distinct:
		[ncd, ind] = unique(ncumsum);
		% interpolate:
		quantiles(:, i) = interp1(ncd, ind/ind(end), quantiles_p);
	end
	
	if p.Results.density
		image_matrix = flipud(results.density');
		if p.Results.normalize_density_plot
			image_matrix = log(1+image_matrix);
			maplength = length(p.Results.density_colormap);
			image_matrix = maplength*(image_matrix/5);
			image(...
				[1 n], ...
				[1 0], ...
				image_matrix);
		else
			imagesc(...
				p.Results.X, ...
				[1 0], ...
				image_matrix);
		end
		
		set(gca, 'YDir', 'normal');
		colormap(p.Results.density_colormap);
		if ~strcmp(p.Results.density_colorbar_location, 'off')
			colorbar('location', p.Results.density_colorbar_location)
		end
		xlim(p.Results.X([1,end]));
    if p.Results.autoclip_density
			if p.Results.outer_quantiles
				ylim([min(quantiles(1, :)) max(quantiles(5, :))]);
			else
				ylim([...
					min((quantiles(1, :) + quantiles(2, :))/2) ...
					max((quantiles(4, :) + quantiles(5, :))/2)]);
			end
		end
		
		hold on
	end
	
	% plot quantiles:
	l_med = plot(p.Results.X,quantiles(3, :), p.Results.median_plotprobs);
	legend_strings{end+1} = 'median';
	legend_handles(end+1) = l_med;
	
	hold on
	if p.Results.quartiles
		l_quart = plot(p.Results.X,quantiles(2, :), p.Results.quartiles_plotprobs);
		legend_strings{end+1} = 'quantiles .25/.75';
		legend_handles(end+1) = l_quart;
		plot(p.Results.X,quantiles(4, :), p.Results.quartiles_plotprobs);
	end
	if p.Results.outer_quantiles
		l_outer = plot(p.Results.X,quantiles(1, :), p.Results.outer_quantiles_plotprobs);
		legend_strings{end+1} = 'quantiles .025/.975';
		legend_handles(end+1) = l_outer;
		plot(p.Results.X,quantiles(5, :), p.Results.outer_quantiles_plotprobs);
	end
end

% plot moments:
if ~isfield(results, 'density') ||  p.Results.add_moments
	if ~isfield(results, 'moments')
		error('results structure does not contain moments')
	end
	% plot expected value (1st moment):
	m1 = results.moments(1, :);
	l_exp = plot(p.Results.X,m1, p.Results.moment1_plotprops);
	legend_strings{end+1} = 'expected val.';
	legend_handles(end+1) = l_exp;
	hold on
	
	if p.Results.std
		if length(results.moments(:,1)) < 2
			error('results structure does not contain 2nd moment')
		end
		variance = results.moments(2, :) - results.moments(1, :).^2;
		stdev = sqrt(variance);
		std_top = m1 + stdev;
		std_bot = m1 - stdev;
		l_std = plot(p.Results.X,std_top, p.Results.moment2_plotprops);
		legend_strings{end+1} = 'standard dev.';
		legend_handles(end+1) = l_std;
		plot(p.Results.X,std_bot, p.Results.moment2_plotprops);
	end
end

ax1 = gca;
if strcmp(p.Results.XAxisLocation, 'top')
	set(ax1, 'XAxisLocation', 'top');
end


if p.Results.bprob
	if ~isfield(results, 'bprob')
		error('results struct does not contain bprob')
	end
	hold on
	col = p.Results.bprob_color;
	if strcmp(p.Results.XAxisLocation, 'top')
		x2loc = 'bottom';
	else
		x2loc = 'top';
	end
	ax2 = axes('Position',get(ax1,'Position'),...
		'XAxisLocation', x2loc,...
		'xticklabel',[], ...
		'YAxisLocation','right',...
		'Color','none',...
		'XColor', col, 'YColor', col);
	bprob = results.bprob;
	bprob(1) = nan;
	bprob(end) = nan;
	l_bprob = line(1:length(bprob), bprob, 'Color', col, 'Parent', ax2);
	legend_strings{end+1} = 'break prob.';
	legend_handles(end+1) = l_bprob;
	xlim(p.Results.X([1,end]));
	ylimits = get(ax1,'YLim');
	yinc = (ylimits(2)-ylimits(1))/5;
	set(ax1,'YTick',[ylimits(1):yinc:ylimits(2)]);
	ylimits = get(ax2,'YLim');
	yinc = (ylimits(2)-ylimits(1))/5;
	set(ax2,'YTick',[ylimits(1):yinc:ylimits(2)]);
end

grid on

if p.Results.show_legend
	legend(legend_handles, ...
		legend_strings{:},  ...
		'Location', p.Results.legend_location);
end

if ~isempty(p.Results.YLim) ylim(p.Results.YLim); end

hold off

