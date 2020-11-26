# itermae

This is tool that parses FASTQ format reads using patterns. 
Specifically, it uses fuzzy regular expressions, so patterns that allow some
degeneracy and using the sequence, not just position, to parse reads.
Then it rebuilds SAM, FASTQ, or FASTA file streams for piping into other tools
or into other files.

It is pretty much just a wrapper to apply fuzzy regex from the 
[`regex`](https://pypi.org/project/regex/)
to sequences in 
[`Biopython`](https://pypi.org/project/biopython/) 
format. That's pretty much it, but it's designed
to be a flexible command line interface to that, for easy parallelization.

# Availability, installation, 'installation'

Options:

1. Use pip to install `itermae`, so 

    python3 -m pip install itermae

1. You can clone this repo, and install it locally. Dependencies are in
    `requirements.txt`, so 
    `python3 -m pip install -r requirements.txt` will install those.
    But if you're not using pip anyways, then you... do you.

1. You can use [Singularity](https://syslab.org) to pull and run a 
    [Singularity image of itermae.py](https://singularity-hub.org/collections/4537), 
    where everything is already installed.
    This is the recommended usage. This image is built with a few other tools,
    like gawk, perl, and parallel, to make command line munging easier.

# Usage

`itermae` is envisioned to be used in a pipe-line where you just got your
FASTQ reads back, and you want to parse them. You can use `zcat` to feed
small chunks into the tool, develop operations that match, filter, and extract
the right groups to assemble the output you want. Then you wrap it it up behind
`parallel` and feed the whole FASTQ file via `zcat` in on standard in.
This parallelizes with a small memory footprint (tune the chunk size), then
you write it out to disk (or stream into another tool?).

Do one thing well, right?

See the jupyter notebook in `demo/`, and the HTML produced from that in that
same folder. That should have some examples and ideas for how to use it.

I believe I'm the only one using this tool, so let me know if you ever try it.
I'd love to hear about it, and would be very eager to help you use it and
try to adapt it to work to your purposes. 

Oh, and this is for BASH shells on Linux/Unix boxes ! I have no idea how
OSX/windows stuff works. Are you unfamiliar with this? If you're at a 
university, ask your librarian. If you're not, look it up online or use the
lessons at Software Carpentries. Or tweet at me about it...

# Caution!

The output group formation and filtering is just using `eval`. This gives
flexibility, but is nowhere near remotely thinking that it would be anywhere
near anything like secure. So this is for use at the command line on your
computer, not web-facing or anything of the sort. Be responsible.
