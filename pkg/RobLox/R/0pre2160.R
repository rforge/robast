## due to a change to .C in 2.16.0
setHook(packageEvent("RobLox", "onLoad"),
        function() RobLox::.setLMfunctions())

.setLMfunctions <- function(){
  if(getRversion() < "2.16.0"){
#    RL <- asNamespace("RobLox")
#    assign(".getA1.locsc", RobLox:::.getA1.locsc.old, envir = RL)
#    assign(".getA2.locsc", RobLox:::.getA2.locsc.old, envir = RL)
#    assign(".getA.loc", RobLox:::.getA.loc.old, envir = RL)
    .getA1.locsc <- .getA1.locsc.old
    assign(".getA2.locsc", RobLox:::.getA2.locsc.old, envir = RL)
    assign(".getA.loc", RobLox:::.getA.loc.old, envir = RL)
    .getA.sc <- .getA.sc.old
    .geta.locsc <- .geta.locsc.old
    .geta.sc <- .geta.sc.old
    .getb.loc <- .getb.loc.old
    .getb.locsc <- .getb.locsc.old
    .getb.sc <- .getb.sc.old
  }
}