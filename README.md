# itermae

Command-line utility to recognize patterns in input sequences and generate 
outputs from groups recognized. Basically, a utility for applying fuzzy regular
expression operations to (primarily) DNA sequence for purposes of DNA 
barcode/tag/UMI parsing, sequence and quality -based filtering, 
and general output re-arrangment.

See the [concept here](https://darachm.gitlab.io/itermae/concept.html).

![itermae diagram](https://darachm.gitlab.io/itermae/_images/parse_diagram_1.svg)

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

1. You can use [Singularity](https://syslab.org) to pull and run a 
    [Singularity image of itermae.py](https://singularity-hub.org/collections/4537), 
    where everything is already installed.
    This is the recommended usage. 

    This image is built with a few other tools,
    like g/mawk, perl, and parallel, to make command line munging easier.

# Usage

`itermae` is envisioned to be used in a pipe-line where you just got your
DNA sequencing FASTQ reads back, and you want to parse them. 
The recommended interface is the YAML config file, as demonstrated
in [the tutorial](https://darachm.gitlab.io/itermae/usage/tutorial.html)
and detailed again in the 
[configuration details](https://darachm.gitlab.io/itermae/usage/config.html).
You can also use a command-line argument interface as detailed more
[in the examples](https://darachm.gitlab.io/itermae/usage/examples.html).

I recommend you test this on small batches of data,
then stick it behind GNU `parallel` and feed the whole FASTQ file via 
`zcat` in on standard input.
This parallelizes with a small memory footprint, then
you write it out to disk (or stream into another tool).

# Thanks

Again, the tool is built upon on the excellent work of 

- [`regex`](https://pypi.org/project/regex/)
- [`Biopython`](https://pypi.org/project/biopython/)
- [`parallel`](https://www.gnu.org/software/parallel/)

# Development, helping

Any issues or advice are welcome as an 
[issue on the gitlab repo](https://gitlab.com/darachm/itermae/-/issues).
Complaints are especially welcome.

For development, see the 
[documentation as rendered from docstrings](https://darachm.gitlab.io/itermae/package.html).

A set of tests is written up with `pytest` module, and can be run from inside
the cloned repo with `make test`.
See `make help` for more options, such as building, installing, and uploading.

There's also a bash script with some longer runs in 
`profiling_tests`, these generate longer runs for profiling purposes
with `cProfile` and `snakeviz`.
But is out of date. Todo is to re-configure and retest that for speed.
