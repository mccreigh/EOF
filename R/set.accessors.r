######################################################################
## a generic accessor function generator

## james mccreight 2 gmail > com
## gpl, no warranties of any kind

## this is meant to operate for just one class/object at a time but defining the accessor function
## for an arbitraty number of slots of that class/object.
## the convention followed here is, for a slot SLOT, get method is a.SLOT and set method is a.SLOT<- .
## the a prefix is customizible but right now is required to be the saem for get and set.

## arguments:
## obj.class is a charatcter string of the class of the object with slots, as in the signature.
## prefix='a.'
## obj.dum is character string for a dummy name of the object with the slots.

## value:
## returns a data fram on the existence of both set and get methods for each slot.

set.accessors <- function( obj.class, prefix='a.', obj.dum='x' ) {

  slots=names(getSlots(obj.class))
  s.a <- as.pairlist(paste(prefix, slots, sep=''))
  names(s.a) <- slots
  
  for (pp in 1:length(s.a)) {
    slot <- names(s.a)[pp]
    a.fun <-s.a[[pp]]
    a.fun.set <- paste(a.fun,"<-",sep='')
    
    ## "get" method
    ## create a generic function
    if (!isGeneric(a.fun)) {  
      fun <- if (is.function(a.fun)) {
        eval(parse(text=slot)) 
      } else {
        eval(parse(text=paste('function(',obj.dum,') standardGeneric(\"',a.fun,'\")',sep='')))
      }
      setGeneric(a.fun, fun)
    }
    
    setMethod(a.fun, obj.class,
              eval(parse(text=paste('function(',obj.dum,') ',obj.dum,'@',slot, sep=''))) )

    ## "set" method
    if (!isGeneric(a.fun.set)) {  
      fun <- if (is.function(a.fun.set)) {
        eval(parse(text=slot)) 
      } else {
        eval(parse(text=paste('function(',obj.dum,', value) standardGeneric(\"',a.fun.set,'\")',sep='')))
      }
      setGeneric(a.fun.set, fun)
    }

    setReplaceMethod(a.fun, obj.class,
                     eval(parse(text=paste('function(',obj.dum,', value)',
                                  ' { ',obj.dum,'@',slot,' <- value ; ',obj.dum,' }',sep='')))
                     )

  }

  all.names <- paste(rep(as.character(s.a),each=2), c('','<-'), sep='')
  return( data.frame( names=all.names, exists=exists( all.names ) ) )
  
}
