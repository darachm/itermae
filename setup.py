import setuptools

setuptools.setup(
    name='itermae',
    version='0.4.0',
    description='desc',
    url='http://gitlab.com/darachm/itermae',
    author='Darach Miller',
    author_email='darachm@stanford.edu',
    license='GPL',
    packages=['itermae'],
    install_requires=[
        'regex',
        'biopython'
        ],
    scripts=['bin/itermae'],
    zip_safe=False
    )
