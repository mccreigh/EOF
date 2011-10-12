n=34
m=20

xnum=2
ynum=0
omega=4

u=seq(0,2*pi,length.out=m)
v=seq(0,2*pi,length.out=n)

y=u %o% rep(1,length(u))
x=t(y)

a=array(data=NA, dim=c(m,m,n))
for (j in 1:n) 
  a[,,j] = 10*sin(ynum*y+ xnum*x - omega*v[j]) + 10*rnorm(length(u));

z=array(data=a, dim=c(m^2,n))

src("ceof.r")
c <- ceof( z )

plot.ceof <- function(x, y, time, c, mode=1 ) {

  vplayout <- function(x,y) viewport(layout.pos.row=x,layout.pos.col=y)
                                     
  ## space
  ceof.space.frame <-
    data.frame( x=as.vector(x), y=as.vector(y),
               real= as.vector(t( Re(c$EOF[,mode]))),
               imaginary= as.vector(t( Im(c$EOF[,mode]))),
               amp= as.vector(t(Mod(c$EOF[,mode]))),
               phase.x= as.vector(t(Re(c$EOF[,mode]))),
               phase.y= as.vector(t(Im(c$EOF[,mode])))
               )
  ## time
  ceof.time.frame <-
    melt(data.frame( time=time,
               real=as.vector( Re( c$proj[,mode] )),
               imaginary=as.vector( Im( c$proj[,mode] ))), id.vars='time')
  
  r <- ggplot(data=ceof.space.frame,
              aes(x=x, y=y, fill=real, z=real )) +
                geom_tile() + stat_contour() 

  i <- ggplot(data=ceof.space.frame,
              aes(x=x, y=y, fill=imaginary, z=imaginary )) +
                geom_tile() + stat_contour() 

  pa <- ggplot(data=ceof.space.frame,
         aes(x=x, y=y, fill=amp, z=amp )) +
         geom_tile() + stat_contour() +
         geom_segment( aes( xend=x+3*phase.x, yend=y+3*phase.y ),
                      arrow=arrow(length=unit(.1,'cm')), color='white' ) 

  t <- ggplot( ceof.time.frame, aes( x=time, y=value, colour=variable) ) +
    geom_line()

  quartz()
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,2)))
  print(r, vp=vplayout(1,1))
  print(i, vp=vplayout(1,2))
  print(pa, vp=vplayout(2,1))
  print(t, vp=vplayout(2,2))
  
}

quartz()
plot.ceof(x, y, 1:n, c , mode=1 )
