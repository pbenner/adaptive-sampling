# Copyright (C) 2011, 2012 Tobias Elze
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#' Plot the density of the binning posterior.
#'
#' @export

plot.binning.density <- function(
	results,
	plot.density = TRUE,
	densitycolor = gray((256:0)/256),
	autoclip.density = TRUE,
	normalize.density.plot = FALSE,
	median.color = 'red',
	median.lty = 'solid',
	plot.quartiles = TRUE,
	quartiles.color = 'red',
	quartiles.lty = 'dashed',
	plot.outer.quantiles = TRUE,
	outer.quantiles.color = 'red',
	outer.quantiles.lty = 'dotted',
	plot.moments = FALSE,
	moment1.color = 'red',
	moment1.lty = 'solid',
	moment2.color = 'red',
	moment2.lty = 'dashed',
	plot.std = TRUE,
	plot.bprob = FALSE,
	bprob.color = 'darkgreen',
	bprob.lty = 'solid',
	bprob.label = "break probabilities",
	show.legend = TRUE,
	legend.location = 'bottomright',	# given either as string or as list or vector of 2 values, e.g. list(0, 1)
	show.grid = TRUE,
	xaxis.location = 'bottom',
	x.labels = NULL,
	xlab = "s",
	ylab=expression(p["x,s"]),
	...
)
{
	legendentries = list(legend='median', col = median.color, lty = list(median.lty))
	addToLegend <- function(l, s, col, lty)
	{
		l$legend = c(l$legend, s)
		l$col = c(l$col, col)
		l$lty[[length(l$lty)+1]] = lty
		l
	}
	isfield <- function(x) { x %in% names(results) }
	plotActive = FALSE
	if(plot.bprob) 
	{
		margins = par("mar")
		if(margins[4] < margins[2]) margins[4] <- margins[2]
		op <- par(mar=margins)
	}
	
	if(isfield('density'))
	{
		plotActive = TRUE
		n <- dim(results$density)[1]
		quantiles.p <- c(.025, .25, .50, .75, .975)
		quantiles <- matrix(numeric(5*n), ncol=n)
		for(i in 1:n)
		{
			v = results$density[i,]
			# normalized cumsum:
			ncumsum = cumsum(v)/sum(v)
			# make data distinct:
			ncd = unique(ncumsum)
			ind = (1:length(ncumsum))[!duplicated(ncumsum)]
			# interpolate:
			interp <- approxfun(ncd, ind/tail(ind,1))
			quantiles[, i] = sapply(quantiles.p, interp)
		}
		
		if(plot.density)
		{
			imagematrix = results$density
			if(normalize.density.plot)
			{
				imagematrix = log(1+imagematrix)
				maplength = length(densitycolor)
				imagematrix = maplength*(imagematrix/5);
			}
			if(is.null(x.labels)) { x.labels = 1:n }
			xaxislabel = ifelse(xaxis.location == 'bottom', xlab, "")
			ylimits = if(autoclip.density) 
				ifelse(
					rep(plot.outer.quantiles, 2), 
					c(max(0, min(quantiles[1,])-0.01), min(1, max(quantiles[5,])+0.01)), 
					c(min((quantiles[1,] + quantiles[2,])/2), max((quantiles[4,] + quantiles[5,])/2)))
				else {NULL}
			image(
				x = x.labels,
				z = imagematrix, 
				col=densitycolor, 
				ylim = ylimits,
				xaxt='n',
				xlab = xaxislabel,
				ylab=ylab,
				...)
			axis(ifelse(xaxis.location == 'bottom', 1, 3), xaxp = par("xaxp"))
			box(bty=ifelse(plot.bprob, 'c', 'o'))

			# plot quantiles:
			lines(x.labels, quantiles[3,], col = median.color, lty = median.lty)
			if(plot.quartiles)
			{
				legendentries = addToLegend(legendentries, '25%/75%', col = quartiles.color, lty = quartiles.lty)
				lines(x.labels, quantiles[2,], col = quartiles.color, lty = quartiles.lty)
				lines(x.labels, quantiles[4,], col = quartiles.color, lty = quartiles.lty)
			}
			if(plot.outer.quantiles)
			{
				legendentries = addToLegend(legendentries, '2.5%/97.5%', col = outer.quantiles.color, lty = outer.quantiles.lty)
				lines(x.labels, quantiles[1,], col = outer.quantiles.color, lty = outer.quantiles.lty)
				lines(x.labels, quantiles[5,], col = outer.quantiles.color, lty = outer.quantiles.lty)
			}
		}
	}
	# plot moments:
	if(!isfield('density') || plot.moments)
	{
		if(!isfield('moments')) { warning("could not find moments", call.=F) }
		else
		{
			m1 = results$moments[1, ]
			if(is.null(x.labels)) { x.labels = 1:length(m1) }
			if(plotActive) { lines(x.labels, m1, col = moment1.color, lty = moment1.lty) }
			else { plot(x.labels, m1, col = moment1.color, type = 'l', lty = moment1.lty, xlab=xlab, ...) }
			legendentries = addToLegend(legendentries, '1st moment', col = moment1.color, lty = moment1.lty)
			
			if(plot.std && dim(results$moments)[1]>1)
			{
				variance = results$moments[2, ] - results$moments[1, ]**2;
				stdev = sqrt(variance);
				std.top = m1 + stdev;
				std.bot = m1 - stdev;
				lines(x.labels, std.top, col = moment2.color, lty = moment2.lty);
				lines(x.labels, std.bot, col = moment2.color, lty = moment2.lty);
				legendentries = addToLegend(legendentries, '2nd moment', col = moment2.color, lty = moment2.lty)
			}
		}
	}
	
	# plot bprob:
	if(plot.bprob)
	{
		if(!isfield('bprob')) { warning("could not find bprob", call.=F) }
		else
		{
			bprob = results$bprob
			bprob[1] = NaN; bprob[length(bprob)] = NaN; 
			yticks = par("yaxp")
			# rescale:
			bprob = bprob*(yticks[2] - yticks[1]) + yticks[1]
			labs = sapply(
				seq(0, 1, length.out = yticks[3]+1), 
				function(s) ifelse(round(s, 1) == s, sprintf("%.1f", s), sprintf("%.2f", s)))
			axis(
				4, 
				at = seq(yticks[1], yticks[2], length.out = yticks[3]+1),
				labels=labs,
				col = bprob.color,
				col.ticks = bprob.color,
				col.axis = bprob.color)
			lines(x.labels, bprob, col = bprob.color, lty = bprob.lty)
			legendentries = addToLegend(legendentries, 'break prob.', col = bprob.color, lty = bprob.lty)
			mtext(bprob.label, side=4, line=3, col=bprob.color)
			par(op)
		}
	}
	if(show.grid) grid(col='gray')
	if(show.legend)
	{
		legendloc = if(length(legend.location)<2) list(legend.location, NULL) else as.list(legend.location)
		ltynames = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
		ltys = sapply(legendentries$lty, function(x) ifelse(is.numeric(x), x, which(ltynames==x)))
		legend(x=legendloc[[1]], y=legendloc[[2]], legend=legendentries$legend, col=legendentries$col, lty=ltys, bg="white")
	}
}

plot.binning <- function(
	results, 
	bar.color="darkgreen", 
	xlab = "",
	posterior.label = "posterior",
	plot.bprob=FALSE, 
	...)
{
	op <- par(mfrow=c(2, 1))
	marginright = ifelse(plot.bprob, 4.1, 0.2)
	par(mar=c(0.1, 4.1, 2.1, marginright))
	plot.binning.density(results, plot.bprob=plot.bprob,  xaxis.location = 'top', ...)
	xaxp = par("xaxp")
	par(mar=c(4.1, 4.1, 0.1, marginright))
	bb = barplot(
		results$mpost, 
		col=bar.color, 
		xlab = xlab, 
		ylab=posterior.label,
		names.arg = 1:length(results$mpost),
		xaxs="i",
		yaxs="r")
	axis(1, at=bb, labels=FALSE)
	grid(nx=0, ny=NULL, col="gray")
	barplot(
		results$mpost, 
		col=bar.color, 
		xlab = xlab, 
		ylab=posterior.label,
		names.arg = 1:length(results$mpost),
		xaxs="i",
		yaxs="r",
		add=TRUE)
	box()
	par(op)
}

#' Plot the  binning posterior.
#'
#' @export
#' @S3Method(plot, binning.posterior)

plot.binning.posterior <- function(posterior, ...)
	plot.binning(unclass(posterior), ...)
