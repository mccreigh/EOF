######################################################################
## method: ggCEOF
## takes a CEOF object and returns a list of ggplot objects.
## use ggplot to plot an CEOF obj/class.

if (!isGeneric("ggCEOF")) {  ## creates a generic function
  fun <- if (is.function("ggCEOF")) EOF else function( ceof, ...) standardGeneric("ggCEOF")
  setGeneric("ggCEOF", fun)
}

## plot=FALSE: TRUE results in a canned plot. this argument may be logical or
##             an integer vector containing the modes to be plotted. nmodes
##             is increased to get any modes specified.

setMethod("ggCEOF", c("CEOF"),
          function( ceof,  nmodes=4, worldmap=FALSE,
                   cut.space.re.im=FALSE, cut.space.amplitude=FALSE,
                   plot=FALSE, melt.time=FALSE,
                   phase.mult=3, unwrap.phase=TRUE )
          {

            require(ggplot2)

            plot=as.vector(plot)+0  ## convert logical to numeric           
            
            nmax.ceof=length(a.modes(ceof))
            nmodes=min(nmodes, nmax.ceof)            
            nmodes=max(nmodes, max(plot))           
            
            ## Create the output list
            ggCEOF=list()

            ## #############
            ## variance frame: variance explaned / eigen values.
            variance <- data.frame( modes=a.modes(ceof),
                                   variance.percent=a.variance.frac(ceof) )
            variance <- ggplot( variance ) 
            ggCEOF[[1]] <- variance
            
            ## divide the ggplot objects into the indivudal modes, with
            ## the space, time, and sequence information for each all together.
print(nmodes)
            for (mm in 1:nmodes) {
              
              ## ############################
              ## spatial frame: real, imaginary, amplitude, phase.real, phase.imag
              nspace=length(a.lon(ceof))
              
              space <- data.frame( lon=a.lon(ceof), lat=a.lat(ceof) )
              space$real <- Re(a.CEOF(ceof)[,mm])
              space$imaginary <- Im(a.CEOF(ceof)[,mm])
              space$amplitude <- Mod(a.CEOF(ceof)[,mm])
              space$phase.re <- space$real  ##these are here so they dont get cut.
              space$phase.im <- space$imaginary

              ## deal with spatial cuts
              cut.space.real=cut.space.re.im+0
              cut.space.imaginary=cut.space.re.im+0
              cut.space.amplitude=cut.space.amplitude+0
              variables <- c('real', 'imaginary', 'amplitude')
              for (ff in 1:length(variables)) {                
                cut.values <- get(paste("cut.space.",variables[ff],sep='')) + 0
                if ( cut.values[1]!=0 ) {
                  values <- space[ variables[ff] ][,1]
                  if (length(cut.values)==1) ## if the number of levels is specified
                    breaks <- c(-1)* max(abs(values)) + ((0:cut.values) *diff( c(-1,1)*max(abs(values)))/cut.values)
                  space[ variables[ff] ] <- cut(values, breaks, include.lowest=TRUE)
                }
              }

              space <- ggplot( space )
              
              ## #############
              ## time data frame: real ts, imaginary ts, amplitude ts, phase ts
              time <-  data.frame( time=a.POSIXct(ceof) )
              time$real <- Re(a.timeseries(ceof)[,mm])
              time$imaginary <- Im(a.timeseries(ceof)[,mm])
              time$amplitude <- Mod(a.timeseries(ceof)[,mm])
              time$phase <- Arg(a.timeseries(ceof)[,mm])
              if (unwrap.phase) time$phase=unwrap.phase(time$phase)
              time <- if (plot[1] | melt.time) ggplot( melt(time,id.vars='time') ) else ggplot(time)
              
              ## #############
              ## list for each mode
              assign(paste("m",mm,sep=''),list(space=space, time=time))
              eval(parse(text=paste('ggCEOF[[',mm+1,']] <- get(paste("m",',mm,',sep=""))',sep='')))
              
            }

            ## name the modes in the output list
            names(ggCEOF) <- c('variance', paste("mode.",1:nmodes,sep=''))

            if (plot[1]!=0) {
              ##  the default plotting
              ## show the fraction variance explained separately, in a new window
              dev.new()
              print(ggCEOF$variance + geom_point( aes(x=modes,y=variance.percent) )  +
                    opts(title=paste('Percent variance explained by mode, CEOF of ',
                           a.variable(ceof), ' (', a.corr.covar(ceof),')',sep=''))                    
                    )

              ## a plot for each desired mode
              for (pp in plot) {

                extras <- scale_fill_gradient() #scale_fill_gradientn(colour = topo.colors(7)) 
                extras.cut <-  scale_fill_brewer(palette="Spectral")

                r <- ggCEOF[[pp+1]]$space + geom_tile( aes(x=lon, y=lat, fill=real) ) +
                       stat_contour( aes(x=lon,y=lat,z=real) ) +
                       labs(fill=paste("CEOF\nMode: ",pp,"\nReal\n",gsub(" ","\n",a.variable(ceof)),sep=''))
                i <- ggCEOF[[pp+1]]$space + geom_tile( aes(x=lon, y=lat, fill=imaginary, z=imaginary) ) +
                       stat_contour( aes(x=lon, y=lat, z=imaginary) ) +
                       labs(fill=paste("CEOF\nMode: ",pp,"\nImaginary\n",gsub(" ","\n",a.variable(ceof)),sep=''))
                if (cut.space.re.im[1]==0) { r <- r + extras; i=i + extras } else
                                           {r <- r + extras.cut; i=i + extras.cut}

                pa <- ggCEOF[[pp+1]]$space +
                  geom_tile(aes(x=lon, y=lat, fill=amplitude)) +
                  stat_contour(aes(x=lon, y=lat, z=amplitude )) +
                  geom_segment( aes_string( x="lon", y="lat",
                                            xend=paste("lon+",phase.mult,"*phase.im"),
                                            yend=paste("lat+",phase.mult,"*phase.re" ),
                                            arrow="arrow(length=unit(.1,'cm'))" ),
                                color='white' )+
                  labs(fill=paste("CEOF\nMode: ",pp,"\nAmplitude\n",gsub(" ","\n",a.variable(ceof)),sep=''))
                if (cut.space.amplitude[1]==0) pa <- pa + extras else pa <- pa + extras.cut

                
                t <- ggCEOF[[pp+1]]$time + geom_line( aes(x=time, y=value) ) +
                       facet_wrap( ~variable, ncol=1, scale='free_y')
                
                dev.new()
                vplayout <- function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)
                grid.newpage()
                pushViewport(viewport(layout=grid.layout(2,2)))
                print(r, vp=vplayout(1,1))
                print(i, vp=vplayout(1,2))
                print(pa, vp=vplayout(2,1))
                print(t, vp=vplayout(2,2))                
              
              }

            }

            invisible(ggCEOF)

          }

          )
