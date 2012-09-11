## due to a change to .C in 2.16.0

.getA1.locsc <- if(getRversion() < "2.16.0"){
                  function(v) .getA1.locsc.old(v)
                }else{
                  function(v) .getA1.locsc.new(v)
                }

.getA2.locsc <- if(getRversion() < "2.16.0"){
                  function(v) .getA2.locsc.old(v)
                }else{
                  function(v) .getA2.locsc.new(v)
                }

.geta.locsc <- if(getRversion() < "2.16.0"){
                  function(v) .geta.locsc.old(v)
               }else{
                  function(v) .geta.locsc.new(v)
               }

.getb.locsc <- if(getRversion() < "2.16.0"){
                  function(v) .getb.locsc.old(v)
               }else{
                  function(v) .getb.locsc.new(v)
               }
                
.getA.sc <- if(getRversion() < "2.16.0"){
               function(v) .getA.sc.old(v)
            }else{
               function(v) .getA.sc.new(v)
            }

.geta.sc <- if(getRversion() < "2.16.0"){
               function(v) .geta.sc.old(v)
            }else{
               function(v) .geta.sc.new(v)
            }

.getb.sc <- if(getRversion() < "2.16.0"){
               function(v) .getb.sc.old(v)
            }else{
               function(v) .getb.sc.new(v)
            }
                
.getA.loc <- if(getRversion() < "2.16.0"){
                function(v) .getA.loc.old(v)
             }else{
                function(v) .getA.loc.new(v)
             }

.getb.loc <- if(getRversion() < "2.16.0"){
                function(v) .getb.loc.old(v)
             }else{
                function(v) .getb.loc.new(v)
             }
