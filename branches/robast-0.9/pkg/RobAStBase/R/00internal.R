#------------------------------------------------------------------------------
# .format.perc : for formatting percentages
#------------------------------------------------------------------------------
### code borrowed from non-exported code from confint.default from package stats
.format.perc <- function (probs, digits)
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
    "%")

