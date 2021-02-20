# itermae 0.4.2

Command-line utility to apply a series of fuzzy regular expression operations
to sequences from a variety of formats, then reconstruct a variety of output
formats from the captured groups -- after applying custom filters on 
matched group position, length, sequence, and/or quality statistics.
Reads and makes FASTQ, FASTA, text-file, and SAM (tab-delimited).
Designed to function with sequence piped in from tools like GNU `parallel`
to permit light-weight parallelization.
Matching is handled as strings in 
[`regex`](https://pypi.org/project/regex/),
and [`Biopython`](https://pypi.org/project/biopython/) is used to represent,
slice, and read/output formats.

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
you write it out to disk (or stream into another tool).

**Tutorial** / **demo**  - there's a jupyter notebook in this root directory
(`demos_and_tutorial_itermae.ipynb`) and the rendered output HTML.
That should have some examples and ideas for how to use it.
There's also some longer runs that are launched by a bash script in
`profiling_tests`, these generate longer runs for profiling purposes
with `cProfile` and `snakeviz`.

Designed for use in command-line shells on a \*nix machine.
