
.PHONY: container run-demos clean dist-pkg upload-pypi

demo: demo/demos_and_tutorial_itermae.html

%.html : %.ipynb
	jupyter nbconvert --to=html --ExecutePreprocessor.timeout=-1 \
		--ExecutePreprocessor.allow_errors=True \
		--execute $<

# This is to clean up after the demo
clean: 
	rm demo/failed.fastq || echo ""
	rm demo/out.sam || echo ""
	rm demo/report.csv || echo ""

pkg-files=setup.py bin/itermae itermae/__init__.py
dist-fils=itermae-0.4.0-py3-none-any.whl  itermae-0.4.0.tar.gz

$dist-files: $pkg-files
	python3 setup.py sdist bdist_wheel

upload-pypi: $dist-files
	python3 -m twine upload --repository pypi $^

itermae.singularity: Singularity
	sudo rm -r $@ || echo "already gone"
	sudo singularity build $@ $<

/tmp/itermae-test-simg: Singularity
	sudo rm -r $@ || echo "already gone"
	sudo singularity build --sandbox $@ $<

