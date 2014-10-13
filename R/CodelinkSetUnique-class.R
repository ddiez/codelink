# TODO: Add comment
# 
# Author: diez
###############################################################################


setClass("CodelinkSetUnique", contains="ExpressionSet")
setMethod("initialize", "CodelinkSetUnique",
function(.Object, sd = new("matrix"), ...)
{
	callNextMethod(.Object,	sd = sd, ...)
})

setValidity("CodelinkSetUnique", function(object)
{
	assayDataValidMembers(assayData(object), c("sd"))
})

