# Demo 01
# 
# This is just a real simple example to start.
# 
# For maximum transparency, in these demos I contrive a fake FASTQ file as
# a string, delimited by the usage of '<<EOF' and 'EOF' as seen below, and feed
# that either directly into the program (as here) or into a temporary file
# (for later demos).

input_fastq=$(cat <<EOF
@NB501157:100:H5J5LBGX2:1:11101:10000:10043 1:N:0:
TTCACGTCCTCGAGGTCTCTTCAGTCGTAGCAGTTCGATGCGTACGCTACAGGTCGACGGTAAGAGAGGGATGTG
+
AAAAAEEEA/<EAEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEE
@NB501157:100:H5J5LBGX2:1:11101:10000:10138 1:N:0:
GCTTCGTCCTCGAGGTCTCTTGGGCAGACACAACGCTACACGTACGCTGCAGGTCGAGGGCACGCGAGAGATGTG
+
AAAAAEEEE/AEAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
EOF
)

echo "${input_fastq}" | 
    ./itermae.py -o "input > (?P<sample>[ATCGN]{5})(?P<after>[ATCGN]{10})" \
        -oid "input.id+'_'+sample.seq" -oseq "after"

# This should output two SAM-format entries, putting the first five bases in the
# read ID and saving the sequence as the next 10 bases after those.

# The above command could be equivalently written with the long-form of the
# flags as --output-id and --output-seq

#		-o "Sample:  input   > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
#		-o "Strain:  rest	 > (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})CGTACGCTGCAGGTCGAC" \
#		-o "UMITail: rest    > (?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}" \
#		-o "UMI:     UMItail > (GAC(?P<umi1>[ATCG])G(?<umi2>[ATCG])A(?<umi3>[ATCG])G(?<umi4>[ATCG])G(?<umi5>[ATCG])G(?<umi6>[ATCG])G){e<=2}" \
#		--output-seq "sample+spacer+strain" \
#		--output-id "input.id+'_umi='+umi1.seq+umi2.seq+umi3.seq+ \
#			umi4.seq+umi5.seq+umi6.seq+'_sample='+sample.seq" \
#		--filter "sample_length == 5 and rest_start >= 16" \
#		--output-id "input.id+'_umi='+umi1.seq+umi2.seq+umi3.seq+ \
#			umi4.seq+umi5.seq+umi6.seq+'_sample='+sample.seq" \
#		--filter "sample_length == 5 and rest_start >= 16" \
#		--job-polling-delay 0.1 \
#		--verbose #-m 
