
## to generate a discrete/factor variable by specifying the step size.
## gives an odd number of levels by default, centered on 0.
## can also return an even number of steps, with zero "between" levels (not actually)
## min(levels)=max(levels)=max(abs(data))

center.discrete <- function( cont, step=NA, odd=TRUE, ...  ) {
  max.abs=max(abs(cont),na.rm=TRUE)
  breaks2 <- ceiling((max.abs-(step/2))/step)
  breaks <- if (odd)  ((-1*breaks2):(breaks2 + 1)-.5)*step else ((-1*breaks2):(breaks2))*step 
  out <- cut( cont, breaks=breaks, ...)
  out <- factor( out, levels=rev(levels(out)))
}
