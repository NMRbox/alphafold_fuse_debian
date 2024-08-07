.DEFAULT_GOAL := build 
.PHONY: build clean wheel

#virtual environments are not relocatable, so we must create it where we want it installed
IDIR := /opt/nmrbox/alphafold_fuse

$(IDIR):
	#create virtual environment
	python3.12 -m venv $(IDIR)

dist/alphafold_fuse*whl: | $(IDIR)
	$(IDIR)/bin/pip install -U pip build
	$(IDIR)/bin/python -m build .

# make target that can be invoked from command line, for testing
wheel: dist/alphafold_fuse*whl

$(IDIR)/bin/alphafold_fuse:  dist/alphafold_fuse*whl
	$(IDIR)/bin/pip install --force-reinstall dist/alphafold_fuse-*.whl

# copy virtual environment into local directory so dh_install can find it
alphafold_fuse:  $(IDIR)/bin/alphafold_fuse
	cp -r $(IDIR) alphafold_fuse 

alphafold_fuse/nmrbox-alphafuse.service:  | alphafold_fuse 
	cp nmrbox-alphafuse.service alphafold_fuse/

alphafold_fuse/update_updatedb.py:  alphafold_fuse/nmrbox-alphafuse.service
	cp update_updatedb.py  alphafold_fuse/

alphafold_fuse/README.md: README.md alphafold_fuse/update_updatedb.py
	cp README.md alphafold_fuse

build: alphafold_fuse/README.md
	
clean:
	rm -fr $(IDIR)  alphafold_fuse dist
	
