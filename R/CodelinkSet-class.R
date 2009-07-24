# CodelinkRawSet-class.
#setClass("CodelinkSet", contains="eSet")
setClass("CodelinkSet", contains="ExpressionSet")

setValidity("CodelinkSet", function(object) {
	#assayDataValidMembers(assayData(object), c("intensity", "background", 
	assayDataValidMembers(assayData(object), c("background", "flag", "snr"))
})
