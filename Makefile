# For help, do `make help` ( idea from SoftwareCarpentry and victoria.dev )

.PHONY: help container run-demos clean dist-pkg dist-files upload-pypi \
	upload-testpypi profiler-runs

help: ## Display help
	@echo 'Commands/rules to run:'
	@grep '\s##\s' Makefile | mawk -F':.*?## ' '{printf "    %-28s%s\n", $$1, $$2}'

demo: demo/demos_and_tutorial_itermae.html ## Typeset the demo jupyter notebook

%.html : %.ipynb
	jupyter nbconvert --to=html --ExecutePreprocessor.timeout=-1 \
		--ExecutePreprocessor.allow_errors=True \
		--execute $<

clean: ## Cleanup the demo outputs 
	rm demo/failed.fastq || echo ""
	rm demo/out.sam || echo ""
	rm demo/report.csv || echo ""

pkg-files=setup.py bin/itermae itermae/__init__.py

test: dist-files ## Run pytest to ... well ... test it
	python3 -m pip install dist/*.whl
	pytest

dist-files: $(pkg-files) ## Make distribution files for pypi
	rm dist/* || echo""
	python3 setup.py sdist bdist_wheel

upload-pypi: dist-files ## Upload distribution files for pypi
	python3 -m twine upload --repository pypi dist/*

upload-testpypi: dist-files ## Upload distribution files for TEST pypi
	python3 -m twine upload --repository testpypi dist/*

profiler-runs: ## Run profiler experiments to look for performance tweaks with snakeviz
	bash profiling_tests/profiler_runs.sh

itermae.singularity: Singularity* ## Build Singularity container from recipe
	sudo rm -r $@ || echo "already gone"
	sudo singularity build $@ $<

itermae-test.singularity: Singularity ## Build TEST sandbox Singularity container from recipe
	sudo rm -r $@ || echo "already gone"
	sudo singularity build --sandbox $@ $<

