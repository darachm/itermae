Basic Usage
============

`itermae` is designed to be a command-line tool to parse FASTQ files.
Here are some basic uses of that, and by that I mean so basic as to be 
impractical for normal use. Normal use should involve the parallelization using
GNU `parallel` detailed in the `parallelization`_ section.

.. _parallelization: parallel


# Demo 01
# 
# This is just a real simple example to start.
# 
# For maximum transparency, in these demos I contrive a fake FASTQ file as
# a string, delimited by the usage of '<<EOF' and 'EOF' as seen below, and feed
# that either directly into the program (as here) or into a temporary file
# (for later demos). Note that for 'echo'ing this variable, the quotes around
# it are required to not munge the newlines.

input_fastq=$(cat <<EOF
@NB501157:100:H5J5LBGX2:1:11101:10000:10043 1:N:0:
TTCACGTCCACGAGGTCTCTTCAGTCGTAGCAGTTCGCGTACGCTACAGGTCGACGGTAAGAGAGGGATGTGATC
+
AAAAAEEEA/<EAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEE
@NB501157:100:H5J5LBGX2:1:11101:10000:19701 1:N:0:
CTACTGTCCACGAGGTCTCTGATGCACTGCGTTCCATGTTCGTACGCTGCAGGTCGACGGAAGGAGCGCGATGTG
+
AAAAAEEEE/AEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
EOF
)

echo "${input_fastq}" | \
    itermae -o "input > (?P<sample>[ATCGN]{5})(?P<after>[ATCGN]{10})" \
        -oid "input.id+'_'+sample.seq" -oseq "after"

# This should output two SAM-format entries, putting the first five bases in the
# read ID and saving the sequence as the next 10 bases after those.
#
# The above command could be equivalently written with the long-form of the
# flags as --output-id and --output-seq
# Demo 02
# 
# This is the same operation as the first demo, but now demonstrating options
# to debug and get verbosity.

input_fastq=$(cat <<EOF
@NB501157:100:H5J5LBGX2:1:11101:10000:10043 1:N:0:
TTCACGTCCACGAGGTCTCTTCAGTCGTAGCAGTTCGCGTACGCTACAGGTCGACGGTAAGAGAGGGATGTGATC
+
AAAAAEEEA/<EAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEE
@NB501157:100:H5J5LBGX2:1:11101:10000:19701 1:N:0:
CTACTGTCCACGAGGTCTCTGATGCACTGCGTTCCATGTTCGTACGCTGCAGGTCGACGGAAGGAGCGCGATGTG
+
AAAAAEEEE/AEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
EOF
)

echo
echo "With level 1 verbosity"
echo "${input_fastq}" | 
    itermae -o "input > (?P<sample>[ATCGN]{5})(?P<after>[ATCGN]{10})" \
        -oid "input.id+'_'+sample.seq" -oseq "after" \
        -v

echo
echo "With level 2 verbosity"
echo "${input_fastq}" | 
    itermae -o "input > (?P<sample>[ATCGN]{5})(?P<after>[ATCGN]{10})" \
        -oid "input.id+'_'+sample.seq" -oseq "after" \
        -v --verbose

echo
echo "With level 3 verbosity"
echo "${input_fastq}" | 
    itermae -o "input > (?P<sample>[ATCGN]{5})(?P<after>[ATCGN]{10})" \
        -oid "input.id+'_'+sample.seq" -oseq "after" \
        -v -v -v

# This should spit out a lot of outputs to standard error, that tell you exactly
# what is going on.
# Demo 02
# 
# Here, we use the same FASTQ input as the previous demo, but we use a more
# complex set of operations. To capture multiple groups. We also use the
# longer flags, for readability.
#
# The capture groups on this one are complex. This amplicon is designed with
# first five bases being sample index, then it's fixed sequence, then it's
# a strain barcode, then fixed sequence, then 6 bases of UMI alternating with
# fixed sequence. Yep. Kinda complex.

input_fastq=$(cat <<EOF
@NB501157:100:H5J5LBGX2:1:11101:10000:10043 1:N:0:
TTCACGTCCACGAGGTCTCTTCAGTCGTAGCAGTTCGCGTACGCTACAGGTCGACGGTAAGAGAGGGATGTGATC
+
AAAAAEEEA/<EAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEE
@NB501157:100:H5J5LBGX2:1:11101:10000:19701 1:N:0:
CTACTGTCCACGAGGTCTCTGATGCACTGCGTTCCATGTTCGTACGCTGCAGGTCGACGGAAGGAGCGCGATGTG
+
AAAAAEEEE/AEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
EOF
)

echo "${input_fastq}" | 
    itermae \
        --operation "input > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
        --operation "rest > (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})CGTACGC" \
        --operation "rest > (?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}" \
        --operation "UMItail > (GAC(?P<umi1>[ATCG])G(?<umi2>[ATCG])A(?<umi3>[ATCG])G(?<umi4>[ATCG])G(?<umi5>[ATCG])G(?<umi6>[ATCG])G){e<=2}" \
        --output-id "input.id+'_umi='+umi1.seq+umi2.seq+umi3.seq+umi4.seq+umi5.seq+umi6.seq+'_sample='+sample.seq" \
        --output-seq "strain" \
        -v

# This should output one SAM-format entry that matches and can form the output
# appropriately.
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
# Demo XX
# 
# Filters
# Filters are evaluated as straight plain python, but where each group is a 
# variable available, and these each have the attributes of 'start', 'end', and
# 'length'.
#
# Here we use a reduced set of operations to keep it sorta simpler.

input_fastq=$(cat <<EOF
@NB501157:100:H5J5LBGX2:1:11101:10000:10043 1:N:0:
TTCACGTCCACGAGGTCTCTTCAGTCGTAGCAGTTCGCGTACGCTACAGGTCGACGGTAAGAGAGGGATGTG
+
AAAAAEEEA/<EAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEE
@NB501157:100:H5J5LBGX2:1:11101:10000:19701 1:N:0:
CTACTGTCCACGAGGTCTCTGATGCACTGCGTTCCATGTTCGTACGCTGCAGGTCGACGGAAGGAGCGCGATGTG
+
AAAAAEEEE/AEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
EOF
)

echo "${input_fastq}" | 
    itermae \
        --operation "input > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
        --operation "rest > (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})(?P<tail>CGTACGC)" \
        --output-id "input.id+'_sample='+sample.seq" \
        --output-seq "strain" \
		--filter "sample.length == 5 and strain.length == 20" 

# Note that the first one is filtered out, not because the 'sample' barcode is 
# not 5 bases at the start, but 'strain' barcode is not exactly 20 bases. It is
# still captured because I specified that 'strain' could capture between 10 and
# 26 bases, but it's not output because the filter statement is not true.


