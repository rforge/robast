bdpPickands <- function(a=2,xi=0.7,orig=TRUE,pr=T, GEVD=T){
## compute BDP of Pickands model
### orig: use original defintion or ours for beta
### pr output as print
### GEVD: TRUE => GEVD FALSE => GPD
   d <- (2*a^xi-1)^(-1/xi)
   if(GEVD){
      p1 <-  exp(-1/a)
      p2 <-  exp(-1/a^2)
      d1 <- exp(-d)
   }else{
      p1 <- 1-1/a
      p2 <- 1-1/a^2
      d1 <- 1-d
   }
   pd <-if (orig) d1 else p1
   if (pr) print(c(a, p1,p2,p2-p1))
   min(p1, 1-p2, p2-pd)
}
if(FALSE){
###GPD:
###### original definition
bdpPickands(a=2,GEVD=F)
##optimal:
ao <- optimize(bdpPickands, interval=c(1,6),
              pr=F, GEVD=F, maximum=TRUE)$max
bdpPickands(a=ao, GEVD=F)
###### new definition
bdpPickands(a=2, orig=F, GEVD=F)
##optimal:
ao <- optimize(bdpPickands, interval=c(1,6),
              orig=F, pr=F, GEVD=F, maximum=TRUE)$max
bdpPickands(a=ao, orig=F, GEVD=F)

##GEVD:
###### original definition
bdpPickands(a=2)
bdpPickands(a=1/log(2))
##optimal:
ao <- optimize(bdpPickands, interval=c(1,6),
              pr=F, maximum=TRUE)$max
bdpPickands(a=ao)
###### new definition
bdpPickands(a=2, orig=F)
bdpPickands(a=1/log(2), orig=F)
##optimal:
ao <- optimize(bdpPickands, interval=c(1,6),
              orig=F, pr=F, maximum=TRUE)$max
bdpPickands(a=ao, orig=F)

}
