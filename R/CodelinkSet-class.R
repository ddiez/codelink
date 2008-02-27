# CodelinkRawSet-class.
#setClass("CodelinkSet", contains="eSet")
setClass("CodelinkSet", contains="ExpressionSet")

#
setMethod("initialize", "CodelinkSet",
function(.Object,
	#phenoData = new("AnnotatedDataFrame"),
	#featureData = new("AnnotatedDataFrame"),
	#experimentData = new("MIAME"),
	#annotation = character(),
	
	#intensity = new("matrix"),
	background = new("matrix"),
	flag = new("matrix"),
	snr = new("matrix"),
	...
) {
	callNextMethod(.Object,
		#phenoData = phenoData,
		#featureData = featureData,
        #experimentData = experimentData,
        #annotation = annotation,

		#intensity = intensity,
		backrgound = background,
		flag = flag,
		snr = snr,
		...
	)
})

setValidity("CodelinkSet", function(object) {
	#assayDataValidMembers(assayData(object), c("intensity", "background", 
	assayDataValidMembers(assayData(object), c("background", "flag", "snr"))
})
