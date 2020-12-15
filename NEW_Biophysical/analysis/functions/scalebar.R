

scalebar <- function(x1, x2, y1, y2, xscale=NULL, yscale=NULL){
	lines(c(x1, x1, x2), c(y2, y1, y1))
	if (!is.null(xscale)) text((x1 + x2)/2, y1, xscale, pos=1)
	if (!is.null(yscale)) text(x1, (y1 + y2)/2, yscale, pos=4)
	# text(125, -69, "2 mV", pos=4)
}


scalebar2 <- function(dx, dy, xscale=NULL, yscale=NULL, pos='bottomleft'){
	lims <- par('usr')
	dxx <- lims[2]-lims[1]
	dyy <- lims[4]-lims[3]
	
	if (pos=='bottomleft'){
		x1 <- lims[1] + dxx/20
		x2 <- x1 + dx
		y1 <- lims[3] + dyy/20 + dy
		y2 <- y1 + dy
	}

	if (pos=='bottomright'){
		x1 <- lims[2] - dxx/20 - dx
		x2 <- x1 + dx
		y1 <- lims[3] + dyy/20 + dy
		y2 <- y1 + dy
	}

	if (pos=='topleft'){
		x1 <- lims[1] + dxx/20
		x2 <- x1 + dx
		y1 <- lims[4] - dyy/20 - dy
		y2 <- y1 + dy
	}

	if (pos=='topright'){
		x1 <- lims[2] - dxx/20 - dx
		x2 <- x1 + dx
		y1 <- lims[4] - dyy/20 - dy
		y2 <- y1 + dy
	}
	
	lines(c(x1, x1, x2), c(y2, y1, y1))
	if (!is.null(xscale)) text((x1 + x2)/2, y1, xscale, pos=1)
	if (!is.null(yscale)) text(x1, (y1 + y2)/2, yscale, pos=4)
	# text(125, -69, "2 mV", pos=4)
}
