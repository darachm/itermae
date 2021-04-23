# For help, do `make help` ( idea from SoftwareCarpentry and victoria.dev )

.PHONY: help container run-demos clean dist-pkg dist-files upload-pypi \
	docs upload-testpypi profiler-runs docker_build docker_push

help: ## Display help
	@echo 'Commands/rules to run:'
	@grep '\s##\s' Makefile | mawk -F':.*?## ' '{printf "    %-28s%s\n", $$1, $$2}'

docs: source ## Build sphinx documentation
	sphinx-build -b html source docs

install: dist-files ## Reinstall from this locally assembled package
	python3 -m pip install dist/*.whl --force-reinstall

test-install: dist-files ## Reinstall from this locally assembled package
	python3 -m pip install dist/*.whl --force-reinstall --no-deps

test: test-install ## Reinstall and run pytest
	pytest

just-test: ## Just run pytest to test the tests
	pytest

pkg-files=setup.py bin/itermae itermae/__init__.py

dist-files: $(pkg-files) ## Make distribution files for pypi
	rm dist/* || echo""
	python3 setup.py sdist bdist_wheel

upload-pypi: dist-files ## Upload distribution files for pypi
	python3 -m twine upload --repository pypi dist/*

upload-testpypi: dist-files ## Upload distribution files for TEST pypi repo
	python3 -m twine upload --repository testpypi dist/*

profiler-runs: ## Run profiler experiments to look for performance tweaks with snakeviz
	bash profiling_tests/profiler_runs.sh

itermae_%.simg : Singularity.% ## Build Singularity container from recipe
	sudo rm -r $@ || echo "already gone"
	sudo singularity build $@ $<

docker_build: ## Build Docker images for itermae:latest and itermae:plus
	docker build . --target itermae -t darachm/itermae
	docker build . --target itermae-plus -t darachm/itermae:plus

docker_push: docker_build ## Push the docker images for itermae:latest, itermae:plus to docker hub
	docker push darachm/itermae
	docker push darachm/itermae:plus

%.simg : Singularity.%
	sudo rm -r $@ || echo "already gone"
	sudo singularity build $@ $<

