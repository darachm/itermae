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
