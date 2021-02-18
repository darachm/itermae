# For help, do `make help`

.PHONY: help container run-demos clean dist-pkg dist-files upload-pypi \
	upload-testpypi profiler-runs

help: ## Display help
	@echo 'Commands/rules to run:'
	@grep '\s##\s' Makefile | mawk -F':.*?## ' '{printf "    %-20s%s", $$1, $$2}'

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

dist-files: $(pkg-files)
	rm dist/* || echo""
	python3 setup.py sdist bdist_wheel

upload-pypi: dist-files
	python3 -m twine upload --repository pypi dist/*

upload-testpypi: dist-files
	python3 -m twine upload --repository testpypi dist/*

# This is for running a few examples on a lot of reads, for profiling
profiler-runs:
	bash profiling_tests/profiler_runs.sh

itermae.singularity: Singularity*
	sudo rm -r $@ || echo "already gone"
	sudo singularity build $@ $<

/tmp/itermae-test-simg: Singularity
	sudo rm -r $@ || echo "already gone"
	sudo singularity build --sandbox $@ $<

