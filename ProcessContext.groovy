def i = 0
def weightBySample = new HashMap<String, Integer>()
def profileBySample = new HashMap<String, Integer>()

def getPropertyVector = { String context ->
	context.collect { it == "N" ? 0.5 : ((it == "A" || it == "T") ? 0 : 1) }
}

def notFirstLine = false

new File("context.txt").splitEachLine("\t") { splitLine ->
	if (notFirstLine) {
		def name = splitLine[0], 
			context = splitLine[1], 
			weight = splitLine[2].toInteger()

		def x = getPropertyVector(context)

		weightBySample.put(name, (weightBySample[name] ?: 0) + weight)

		x.eachWithIndex { it, ind ->
			def signature = "$name\t$ind".toString()
			profileBySample.put(signature, 
				(profileBySample[signature] ?: 0) + it * weight)
		}
	} else {
		notFirstLine = true
	}
}

new File("context.proc.txt").withPrintWriter { pw ->
	pw.println("name\tpos\tvalue\tsum")

	profileBySample.each {
		def name = it.key.split("\t")[0]
		pw.println(it.key + "\t" + it.value + "\t" + weightBySample[name])
	}
}