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
    ./itermae.py \
        --operation "input > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}" \
        --operation "rest > (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})(?P<tail>CGTACGC)" \
        --output-id "input.id+'_sample='+sample.seq" \
        --output-seq "strain" \
		--filter "sample.length == 5 and strain.length == 20" 

# Note that the first one is filtered out, not because the 'sample' barcode is 
# not 5 bases at the start, but 'strain' barcode is not exactly 20 bases. It is
# still captured because I specified that 'strain' could capture between 10 and
# 26 bases, but it's not output because the filter statement is not true.


