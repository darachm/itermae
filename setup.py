import setuptools

with open('README.md','r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='itermae-darachm',
    version='0.4.0',
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
        'regex',
        'biopython',
        ],
    scripts=['bin/itermae'],
    zip_safe=False,
    python_requires='>=3.6'
    )
