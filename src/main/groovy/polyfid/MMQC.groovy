package polyfid

import com.milaboratory.core.alignment.AffineGapAlignmentScoring
import com.milaboratory.core.alignment.Aligner
import com.milaboratory.core.io.sequence.fastq.PairedFastqReader
import com.milaboratory.core.mutations.MutationType
import com.milaboratory.core.sequence.NSequenceWithQuality
import com.milaboratory.core.sequence.NucleotideSequence

def ref = new NucleotideSequence("TAGCGTGAAGACGACAGAACCATGGGATCCATTATCGGCGGCGAATTTACCACCATTGAAAACCAGCCGTGGTTTGCGGCGA" +
        "TTTATCGTCGTCATCGTGGCGGCAGCGTGACCTATGTGTGCGGCGGCAGCCTGATTAGCCCGTGCTGG")


byte maxQual = 40
float minCoverage = 0.95, minSimilarity = 0.90

def fromToPositionQual = new int[maxQual + 1][][][]

for (int i = 0; i <= maxQual; i++) {
    fromToPositionQual[i] = new int[ref.size()][4][4]
}

def fastq1 = args[0], fastq2 = args[1]

def reader = new PairedFastqReader(fastq1, fastq2)

def scoring = AffineGapAlignmentScoring.getNucleotideBLASTScoring()

def addQual = { int qual, int pos, int from, int to ->
    for (int j = 0; j <= qual; j++) {
        fromToPositionQual[j][pos][from][to]++
    }
}

def alignAndUpdate = { NSequenceWithQuality nsq ->
    def aln = Aligner.alignLocalAffine(scoring, ref, nsq.sequence)

    if (aln != null && (aln.sequence2Range.length() / (float) nsq.sequence.size()) >= minCoverage &&
            aln.similarity() >= minSimilarity) {

        def mutations = aln.absoluteMutations

        def mutatedPositions = new HashSet<Integer>()

        (0..<mutations.size()).each { int i ->
            if (mutations.getTypeByIndex(i) == MutationType.Substitution) {
                int pos = mutations.getPositionByIndex(i),
                    from = mutations.getFromAsCodeByIndex(i),
                    to = mutations.getToAsCodeByIndex(i)

                if (to < 4) {
                    int posInRead = aln.convertPosition(pos)

                    if (posInRead >= 0) {
                        byte qual = Math.min(maxQual, nsq.quality.value(posInRead))

                        addQual(qual, pos, from, to)
                    }
                }

                mutatedPositions.add(pos)
            }
        }

        (aln.sequence1Range.lower..<aln.sequence1Range.upper).each { int pos ->
            if (!mutatedPositions.contains(pos)) {
                int posInRead = aln.convertPosition(pos)

                if (posInRead >= 0) {
                    byte from = ref.codeAt(pos)
                    byte qual = Math.min(maxQual, nsq.quality.value(posInRead))

                    addQual(qual, pos, from, from)
                }
            }
        }
    }
}

def read

int n = 0

while ((read = reader.take()) != null) {
    alignAndUpdate(read.r1.data)
    alignAndUpdate(read.r2.data.reverseComplement)

    if (++n % 100000 == 0) {
        println("[" + new Date().toString() + "] Processed $n reads")
    }
}

println("[" + new Date().toString() + "] Finished processing $n reads, writing output")

new File(args[2]).withPrintWriter { pw ->
    pw.println("qual\tpos\tfrom\tto\tcount\tfreq")

    for (int i = 0; i < maxQual; i++) {
        for (int j = 0; j < ref.size(); j++) {
            for (byte k = (byte) 0; k < (byte) 4; k++) {
                for (byte l = (byte) 0; l < (byte) 4; l++) {
                    double count = fromToPositionQual[i][j][k][l],
						   total = fromToPositionQual[i][j][k][k]
					if (total > 0) {
						pw.println(i + "\t" +
								j + "\t" +
								NucleotideSequence.ALPHABET.codeToSymbol(k) + "\t" +
								NucleotideSequence.ALPHABET.codeToSymbol(l) + "\t" +
								count + "\t" +
								(count / (double) total))
					}
                }
            }
        }
    }
}