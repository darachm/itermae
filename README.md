# itermae

Command-line utility to recognize patterns in input sequences and generate 
outputs from groups recognized. Basically, a utility for applying fuzzy regular
expression operations to (primarily) DNA sequence for purposes of DNA 
barcode/tag/UMI parsing, sequence and quality -based filtering, 
and general output re-arrangment.

Reads and makes FASTQ, FASTA, text-file, and SAM (tab-delimited).
Designed to function with sequence piped in from tools like GNU `parallel`
to permit light-weight parallelization.
Matching is handled as strings in 
[`regex`](https://pypi.org/project/regex/),
and [`Biopython`](https://pypi.org/project/biopython/) is used to represent,
slice, and read/output formats.

Designed for use in command-line shells on a \*nix machine.

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
DNA sequencing FASTQ reads back, and you want to parse them. 

You feed small chunks of the file into the tool with match-level
verbosity and record-level reports to develop good patterns. 
These patterns, filtering, and outputs are used to pull out and 
assemble the output you want.

Then you wrap it it up behind
`parallel` and feed the whole FASTQ file via `zcat` in on standard input.
This parallelizes with a small memory footprint (will measure later), then
you write it out to disk (or stream into another tool).

**Tutorial** / **demo**  - there's a jupyter notebook in this root directory
(`demos_and_tutorial_itermae.ipynb`) and the rendered output HTML.
That should have some examples and ideas for how to use it.
There's also some longer runs that are launched by a bash script in
`profiling_tests`, these generate longer runs for profiling purposes
with `cProfile` and `snakeviz`.

