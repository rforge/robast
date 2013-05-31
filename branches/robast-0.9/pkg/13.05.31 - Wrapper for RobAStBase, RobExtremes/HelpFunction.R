 myplot <- function(x,y, ..., withCall = TRUE){
    ###
    ### 1. grab the dots (and probably manipulate it within the wrapper function)
    ###
  	mc <- as.list(match.call(expand.dots = FALSE))[-1]
  	dots <- mc$"..."
  	if(is.null(mc$withCall)) mc$withCall <- TRUE

  	if(missing(x)) stop("Argument 'x' must be given as argument to 'myplot'")
  	if(missing(y)) stop("Argument 'y' must be given as argument to 'myplot'")
    ###
    ## do something to fix the good default arguments
    ###
    ### 2. build up the argument list for the (powerful/fullfledged)
    ### graphics/diagnostics function;
    ### mind not to evaluate the x and (possibly) y args to provide automatic
    ### axis annotation
    ###
    args <- c(list(x=substitute(x),y=substitute(y)),dots,type="l")
	print(args)
	print("###################################################")
    ###
    ### 3. build up the call but grab it and write it into an object
    ###
    cl <- substitute(do.call(plot,args0), list(args0=args))
	print(cl)
	print("###################################################")
  	### manipulate it so that the wrapper do.call is ommitted
    cl0 <- as.list(cl)[-1]
	print(cl0)
	print("###################################################")
  	mycall <- c(cl0[1],unlist(cl0[-1]))
	print(mycall)
	print("###################################################")
  	mycall <- as.call(mycall)
	print(mycall)
	print("###################################################")
  	###
    ### 4. evaluate the call (i.e., produce the graphic)
    ###
    eval(mycall)
    ###
    ### 5. return the call (if withCall==TRUE)
    ###
    if(mc$withCall) print(mycall)

}

x <- 1:20
y <- rnorm(20)
cl <- myplot(x,y,col="red", withCall=TRUE)
cl <- myplot(x,y,col="blue")
cl <- myplot(x,y,col="green", withCall=FALSE)