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
