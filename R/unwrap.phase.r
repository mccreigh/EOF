## Unwrap phase angles.  Algorithm minimizes the incremental phase variation 
## by constraining it to the range [-pi,pi]

unwrap.phase <- function( p,cutoff=pi) {

  np <- length(p)  
  pdiff <- p[2:np]-p[1:(np-1)]
  pdiff.mod <- ( (pdiff+pi) %% 2*pi ) - pi ## shift phase differences to to [-pi,pi)
  pdiff.mod[which(pdiff.mod==-pi & pdiff>0)] <- pi    ## now in [-pi,pi] 
  unwrap <- pdiff.mod - pdiff              
  unwrap[which(abs(pdiff)<cutoff)] <- 0  ## Ignore correction when incr. variation is < CUTOFF
  
  ## Integrate corrections and add to P to produce smoothed phase values
  p[2:np] = p[2:np] + cumsum(unwrap)
  p

}
