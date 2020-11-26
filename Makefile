
.PHONY: container run-demos clean-demo clean

run-demos: demos_and_tutorial_itermae.ipynb
	jupyter nbconvert --to=html --ExecutePreprocessor.timeout=-1 --execute $^

dist_pkg:
	python3 setup.py sdist bdist_wheel

clean: clean-demo
	rm -r tmp || echo ""
	rm -r .nextflow* || echo ""
	rm -r nextflow* || echo ""
	rm -r work || echo ""

clean-demo:
	rm failed.fastq || echo ""
	rm out.sam || echo ""
	rm report.txt || echo ""
	
itermae.singularity: Singularity
	sudo rm -r $@ || echo "already gone"
	sudo singularity build $@ $<

/tmp/itermae-test-simg: Singularity
	sudo rm -r $@ || echo "already gone"
	sudo singularity build --sandbox $@ $<

