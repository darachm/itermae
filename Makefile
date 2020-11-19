
.PHONY: container run-demos clean

demos=$(wildcard demos/*sh)

run-demos: $(demos)
	mkdir -p demos/tmp
	echo $^ | xargs bash
	

dist_pkg:
	python3 setup.py sdist bdist_wheel

#example: tmp/barseq_shortrun_pass.fastq tmp/barseq_longrun_pass.fastq

#container_test: tmp/barseq_shortrun_pass_container.fastq tmp/barseq_longrun_pass_container.fastq

clean: 
	rm -r tmp || echo ""
	rm -r .nextflow* || echo ""
	rm -r nextflow* || echo ""
	rm -r work || echo ""
	
itermae.singularity: Singularity
	sudo rm -r $@ || echo "already gone"
	sudo singularity build $@ $<

/tmp/itermae-test-simg: Singularity
	sudo rm -r $@ || echo "already gone"
	sudo singularity build --sandbox $@ $<

