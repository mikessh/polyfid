import com.milaboratory.core.alignment.AffineGapAlignmentScoring
import com.milaboratory.core.alignment.Aligner
import com.milaboratory.core.io.sequence.fastq.PairedFastqReader
import com.milaboratory.core.mutations.MutationType
import com.milaboratory.core.sequence.NSequenceWithQuality
import com.milaboratory.core.sequence.NucleotideSequence

@Grapes(
        @Grab(group = 'com.milaboratory', module = 'milib', version = '1.5')
)

def ref = new NucleotideSequence("TAGCGTGAAGACGACAGAACCATGGGATCCATTATCGGCGGCGAATTTACCACCATTGAAAACCAGCCGTGGTTTGCGGCGA" +
        "TTTATCGTCGTCATCGTGGCGGCAGCGTGACCTATGTGTGCGGCGGCAGCCTGATTAGCCCGTGCTGG")


byte maxQual = 40
float minCoverage = 0.8, minSimilarity = 0.8

def fromToPositionQual = new int[maxQual + 1][][][]

for (int i = maxQual; i >= 0; i--) {
    def arr = new int[ref.size()][4][4]

    for (int j = i; j <= maxQual; j++) {
        fromToPositionQual[j] = arr // cumulative
    }
}

def fastq1 = args[0], fastq2 = args[1]

def reader = new PairedFastqReader(fastq1, fastq2)

def scoring = AffineGapAlignmentScoring.getNucleotideBLASTScoring()

def alignAndUpdate = { NSequenceWithQuality nsq ->
    def aln = Aligner.alignLocalAffine(scoring, ref, nsq.sequence)

    if (aln != null && (aln.sequence1Range.length() / (float) ref.size()) > minCoverage &&
            aln.similarity() > minSimilarity) {

        def mutations = aln.absoluteMutations

        def mutatedPositions = new HashSet<Integer>()

        (0..<mutations.size()).each { int i ->
            if (mutations.getTypeByIndex(i) == MutationType.Substitution) {
                int pos = mutations.getPositionByIndex(i),
                    from = mutations.getFromAsCodeByIndex(i),
                    to = mutations.getToAsCodeByIndex(i)

                int posInRead = aln.convertPosition(pos)

                if (posInRead >= 0) {
                    byte qual = Math.min(maxQual, nsq.quality.value(posInRead))

                    fromToPositionQual[qual][pos][from][to]++
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

                    fromToPositionQual[qual][pos][from][from]++
                }
            }
        }
    }
}

def read

while ((read = reader.take()) != null) {
    alignAndUpdate(read.r1.data)
    alignAndUpdate(read.r2.data.reverseComplement)
}

new File(args[2]).withPrintWriter { pw ->
    pw.println("qual\tpos\tfrom\tto\tfreq")

    for (int i = 0; i < maxQual; i++) {
        for (int j = 0; j < ref.size(); j++) {
            for (byte k = (byte) 0; k < (byte) 4; k++) {
                for (byte l = (byte) 0; l < (byte) 4; l++) {
                    double count = fromToPositionQual[i][j][k][l]
                    pw.println(i + "\t" +
                            j + "\t" +
                            NucleotideSequence.ALPHABET.codeToSymbol(k) + "\t" +
                            NucleotideSequence.ALPHABET.codeToSymbol(l) + "\t" +
                            count + "\t" +
                            (count / (double) fromToPositionQual[i][j][k][k]))
                }
            }
        }
    }
}