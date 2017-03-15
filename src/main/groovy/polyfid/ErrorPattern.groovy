package main.groovy.polyfid

def referenceFileName = "data/references.txt", outputFileName = "data/kmer_errors.txt", rebuildRef = false
def p1 = "polerr73",
    p2 = "polerr82"
//def p1 = "polerr2016-1",
//    p2 = "polerr2016-2"

def inputFileNamesMap = [
        "data/${p1}.encyclo.variant.caller.txt"        : "Encyclo",
        "data/${p2}.encyclo.variant.caller.txt"        : "Encyclo",
        "data/${p1}.kappa-hf-taq.variant.caller.txt"   : "Kappa HF",
        "data/${p2}.kappa-hf-taq.variant.caller.txt"   : "Kappa HF",
        "data/${p1}.sd-hs.variant.caller.txt"          : "SD-HS",
        "data/${p2}.sd-hs.variant.caller.txt"          : "SD-HS",
        "data/${p1}.snp-detect.variant.caller.txt"     : "SNP-detect",
        "data/${p2}.snp-detect.variant.caller.txt"     : "SNP-detect",
        "data/${p1}.taq-hs.variant.caller.txt"         : "Taq-HS",
        "data/${p2}.taq-hs.variant.caller.txt"         : "Taq-HS",
        "data/${p1}.tersus.variant.caller.txt"         : "Tersus",
        "data/${p2}.tersus.variant.caller.txt"         : "Tersus",
        "data/${p1}.tersus-snp-buff.variant.caller.txt": "Tersus-SNP-buffer",
        "data/${p2}.tersus-snp-buff.variant.caller.txt": "Tersus-SNP-buffer",
        "data/${p1}.truseq.variant.caller.txt"         : "TruSeq",
        "data/${p2}.truseq.variant.caller.txt"         : "TruSeq",
        "data/${p1}.velox.variant.caller.txt"          : "Velox",
        "data/${p2}.velox.variant.caller.txt"          : "Velox"
]

def getProject = { String str ->
    str.split('[/.]')[1]
}

def parseMutation = { String str ->
    def strSplit = str[1..-1].split("[:>]")

    if (str[0] != "S") return null

    [("pos"): strSplit[0].toInteger(), ("from"): strSplit[1], ("to"): strSplit[2]]
}

def referenceMap = new File(referenceFileName).readLines().collectEntries {
    def splt = it.split("\t")
    [(splt[0]): splt[1]]
}

if (rebuildRef) {
    def refMap2 = new HashMap<String, List<String>>()

    inputFileNamesMap.keySet().each { inputFileName ->
        def lineNum = 0
        new File(inputFileName).splitEachLine("\t") { splitLine ->
            if (++lineNum > 2) {
                def ref = splitLine[0], mut = parseMutation(splitLine[1])

                if (mut) {
                    def refSeqLst = refMap2[ref] ?: new ArrayList<String>(("N" * 1000).split("") as List)

                    refSeqLst.add(mut.pos, mut.from)

                    refMap2[ref] = refSeqLst
                }
            }
        }
    }

    refMap2.each { referenceMap[it.key] = it.value.join("") }
}

def revCompMap = [((char) "A"): "T", ((char) "T"): "A", ((char) "G"): "C", ((char) "C"): "G", ((char) "N"): "N"]

def revComp = { String str -> str.reverse().toCharArray().collect { revCompMap[it] }.join("") }

def getKmer = { String refSeq, String from, int pos, int offsetLeft, int offsetRight ->
    assert refSeq[pos] == "N" || refSeq[pos] == from

    def posFrom = pos - offsetLeft, posTo = pos + offsetRight
    def leftN = "", rightN = ""

    if (posFrom < 0) {
        leftN = "N" * (-posFrom)
        posFrom = 0
    }

    if (posTo >= refSeq.length()) {
        rightN = "N" * (posTo - refSeq.length() + 1)
        posTo = refSeq.length() - 1
    }

    leftN + refSeq[posFrom..posTo] + rightN
}

def getKmerBothStrands = { String kmer ->
    def rcKmer = revComp(kmer)

    [kmer, rcKmer].sort().join(",")
}

new File(outputFileName).withPrintWriter { pw ->
    pw.println("project\treference\tname\tmut.pos\tmut.from\tregion.gc\tcount\tcoverage")
    inputFileNamesMap.each { inputFileEntry ->
        def inputFileName = inputFileEntry.key
        println inputFileName
        def lineNum = 0
        new File(inputFileName).splitEachLine("\t") { splitLine ->
            if (++lineNum > 2) {
                def ref = splitLine[0], mut = parseMutation(splitLine[1]),
                    count = splitLine[2], coverage = splitLine[4] // different formats
                if (mut && (count.toDouble() / coverage.toDouble() < 0.25)) {
                    def refSeq = referenceMap[ref]

                    ([[7, 7]]).each { offsets ->
                        def kmer = getKmer(refSeq, mut.from, mut.pos, offsets[0], offsets[1])

                        def kmer2 = kmer[0..5] + "N" + kmer[7..-1]

                        pw.println([getProject(inputFileName),
                                    ref,
                                    inputFileEntry.value,
                                    mut.pos,
                                    mut.from,
                                    kmer2.replaceAll("[ATN]", "").length() / (float) kmer2.replaceAll("[N]", "").length(),
                                    count, coverage
                        ].join("\t"))
                    }
                }
            }
        }
    }
}

/*
new File(outputFileName).withPrintWriter { pw ->
    pw.println("reference\tname\tmut.pos\tmut.from\tkmer\tgood.kmer\tcount\tcoverage")
    inputFileNamesMap.each { inputFileEntry ->
        def inputFileName = inputFileEntry.key
        println inputFileName
        def lineNum = 0
        new File(inputFileName).splitEachLine("\t") { splitLine ->
            if (++lineNum > 2) {
                def ref = splitLine[0], mut = parseMutation(splitLine[1]),
                    count = splitLine[2], coverage = splitLine[4] // different formats
                if (mut && (count.toDouble() / coverage.toDouble() < 0.25)) {
                    def refSeq = referenceMap[ref]

                    ([[1, 1]]).each { offsets ->
                        def kmer = getKmer(refSeq, mut.from, mut.pos, offsets[0], offsets[1]),
                            kmerExt = getKmerBothStrands(kmer)

                        def kmerExt2 = kmer[0] + "N" + kmer[2]

                        pw.println([getProject(inputFileName),
				    ref,
                                    inputFileEntry.value,
                                    mut.pos,
                                    mut.from,
                                    kmerExt2,//kmerExt,
                                    !kmer.contains("N"),
                                    count, coverage//,
                                    //refSeq.replaceAll("[AT]", "").length() / (float) refSeq.length()
                        ].join("\t"))
                    }
                }
            }
        }
    }
}*/
