###########################################################
## function with Control as S4 class   PR 21.5.12
###########################################################


## utilities:

## we use a slightly different semantics from R here:
## if formals are endowed with default values
##    and if they are not matched by exact or partial matching
##    they get set to default __before__ positional matching

## .fix.in.defaults is to perform the default setting
## (before positional matching)

.fix.in.defaults <- function(call.list, fun){

    formals.fun <- formals(fun)
    k <- length(call.list)
    L <- length(formals.fun)
    if("..." %in% names(formals.fun)) L <- L-1
    for(i in 1:L){
        if(!is(formals.fun[[i]],"name")){
           if(!names(formals.fun)[i] %in% names(call.list)&&
              !is.null(formals.fun[[i]])){
              k <- k + 1
              call.list[[k]] <- formals.fun[[i]]
              names(call.list)[k] <- names(formals.fun)[i]
           }
        }
     }
    return(call.list)

}

## as our Control argument is coming as a (NAMED!) list
## and we want to be able to use these names, we attach the
## list; at the end of the call (with on.exit(), so that this
## is done even if an error was thrown), we call .clean.search
## to remove no-longer needed positions in the search path.

.clean.search<-function(what="structure\\(list\\("){

   se <- search()
   se0 <- grep(what,search(), value=TRUE)
   du <- sapply(se0, detach, character.only=TRUE)
   return(invisible(NULL))
}

## our S4 class FctWithCtrl: the function is _not_ a slot but rather
##    the class it self ( "is" instead of "has" )
##    and has additional slot Ctrl, a list
## to be able to have access to Ctrl, one would need something
## like a "self" which is not there in S4 classes
## so instead we build something around it (see below)

setClass("FctWithCtrl", representation(Ctrl="list"), contains="function")

### generating function generates an object of class FctWithCtrl
## (for simplicity of discussion assumed assigned to symbol "fu")
##     takes two arguments f (the function to be extended) and Ctrl, the
##     control list; this achieves
## (a) [almost] arbitrary functions f can be used as "function" (to be enhanced
##              by Ctrl [restriction: no name of the Ctrl list may appear
##              in the formals nor in the actual arguments (passed through "...")
## (b) an arbitrary list as Ctrl argument; no checking so far...
## (c) a call to fu then behaves exactly as a call to argument f, except that
##     in the body of fu one has access to the items of Ctrl
##     -> no clever handling of masking so far

## trick: we write to layers of wrapping functions which are:
##       + an anonymous top layer function doing the management of
##         call, names etc
##       + an inner function f1 which is a variation of f
##         with a modified arg list (the formals of f1 are the union of formals
##             of f and items of Ctrl)
##         with the same return value as f
##         with the body of f1 being the body of f with additional access
##         to the items of Ctrl
###      + the real trick is that Ctrl is an attribute to the generated
##         object and that we can access these arguments from an inner
##         function with sys.function(1)

FctWithCtrl <- function(f, Ctrl){

     g <- new("FctWithCtrl",
              Ctrl = Ctrl,
              function(...){

                 ## do  exact & partial matching by names
                 mc <- match.call(expand=TRUE)

                 dots <- mc$...
                 mc$... <- NULL

                 Ctrl<- attributes(sys.function(1))$Ctrl
                 nmsall <- unique(c(names(Ctrl),names(formals(f))))

                 ## remove items from dots -> needed?
                 if(length(dots)&&!is.null(names(dots))){
                    for(i in 1: length(dots))
                        if(names(dots)[i] %in% nmsall)  dots[[i]] <- NULL
                 }

                 ## remove items from Ctrl if as named arguments in mc
                 if(length(mc)&&!is.null(names(mc))){
                    for(i in 1: length(Ctrl))
                        if(names(Ctrl)[i] %in% names(mc))  Ctrl[[i]] <- NULL
                 }

                 ## combine to new call
                 mc <- as.call(c(as.list(mc),Ctrl,dots))

                 ## after exact and partial matching and before positional
                 ## matching set formals to defaults:
                 mc <- as.call(.fix.in.defaults(as.list(mc),f))

                 ## positional matching
                 mcfit <- as.list(mc[-1])
                 nmcfit <- names(mcfit)
                 unfit <- setdiff(names(formals(f)),nmcfit)
                 for(i in 1:length(nmcfit)){
                     if(nmcfit[i]=="" && length(unfit)>0){
                        names(mcfit)[i] <- if(unfit[1]!="...") unfit[1] else  ""
                        unfit <- unfit[-1]
                     }
                 }

                 ## produce new call from the manipulated list
                 mc0 <- as.list(mc)
                 for(i in 1: length(mcfit))
                    {mc0[[i+1]] <- mcfit[[i]]
                     names(mc0)[i+1] <- names(mcfit)[i]
                 }
                 mc <- as.call(mc0)

                 ## generate variant f1 from f
                 ##
                 f1 <- f
                 ## combine formals
                 formals(f1) <- c(formals(f),Ctrl)
                 ## manipulate body:

                 body(f1) <- substitute({
                         on.exit(.clean.search())
                         attach(ctrl)
                         ..z <- body.f
                         return(..z)
                         }, list(lu=attributes(sys.function(1)),
                                 ctrl=attributes(sys.function(1))$Ctrl,
                                 body.f = body(f))
                         )
                 ## fall this f1 with the manipulated call
                 erg <- do.call(f1, as.list(mc[-1]))
                 return(erg)
                 })
     return(g)
}


### EXAMPLE:

## we want to produce a function f where all access to Ctrl arguments
##   is demonstrated: (attention
##

f.example <- function(x, z=2, ...){
        ## show the manipulated matched call with items from Ctrl attached
        mc0 <- match.call()
        ## show that we have access to items of Ctrl
        print(c(a=a,A=A))
        ## show the ... argument:
        if(length(list(...)))
           cat(" ...", paste(names(list(...)),collapse=", "), "\n",
               "...", paste(..., collapse=", "),"\n");
        ##show how formals are matched:
        print(c(x=x,z=z));
        ## return a value
        sin(x)}

### attention does not work so far as a and A are not yet bound:
# > f.example(2,3)
# Error in print(c(a = a, A = A)) : object 'a' not found

## fu as S4 object
fu <- FctWithCtrl( f=f.example, Ctrl=list(a=2,A=5))

## calls to fu
fu(2,3)
fu(x=22)
fu(AA=3,3,.a=3)
#> ## calls to fu
#> fu(2,3)
#a A
#2 5
# ...
# ... 3
#x z
#2 2
#[1] 0.9092974
#> fu(x=22)
#a A
#2 5
# x  z
#22  2
#[1] -0.008851309
#> fu(AA=3,3,.a=3)
#a A
#2 5
# ... AA, .a
# ... 3 3
#x z
#3 2
#[1] 0.14112
#>

##manipulate Ctrl
fu@Ctrl <- list(A=4,a=3)
fu(AA=3,3,.a=3)

#> fu@Ctrl <- list(A=4,a=3)
#> fu(AA=3,3,.a=3)
#a A
#3 4
# ... AA, .a
# ... 3 3
#x z
#3 2
#[1] 0.14112
