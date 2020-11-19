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
