# CodelinkSet-class.
setClass("CodelinkSet", contains="ExpressionSet")

setValidity("CodelinkSet", function(object) {
	assayDataValidMembers(assayData(object), c("background", "flag", "snr"))
})
