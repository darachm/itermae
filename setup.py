import setuptools

with open('README.md','r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='itermae',
    version='0.6.0.1',
    author='Darach Miller',
    description='Commandline tool for parsing NGS reads by multiple fuzzy '+
        'regex operations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='http://gitlab.com/darachm/itermae',
    author_email='darachm@stanford.edu',
    license='BSD 2-clause',
    packages=setuptools.find_packages(),
    install_requires=[
        'pyyaml',
        'regex',
        'biopython',
        ],
    scripts=['bin/itermae'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    zip_safe=False,
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        ],
    keywords='fastq regex fuzzy amplicon parser barcode extractor extracter'
    )
