.DEFAULT_GOAL := build 
.PHONY: build  clean wheel

#virtual environments are not relocatable, so we must create it where we want it installed
IDIR := /opt/nmrbox/alphafold_fuse

$(IDIR):
	#create virtual environment
	python3.8 -m venv $(IDIR) 

dist/pipmanager*whl: | $(IDIR)
	$(IDIR)/bin/pip install -U pip build
	$(IDIR)/bin/python -m build .

# make target that can be invoked from command line, for testing
wheel: dist/pipmanager*whl

$(IDIR)/bin/pipmanager:  dist/pipmanager*whl
	$(IDIR)/bin/pip install --force-reinstall dist/pipmanager-*.whl

# we can't (or don't know how to) override rules because we have two debian/*install files,
# so copy virtual environment into local directory should dh_install can find it
pipmanager:  $(IDIR)/bin/pipmanager 
	cp -r $(IDIR) pipmanager

build: pipmanager 
	
clean:
	rm -fr $(IDIR)  pipmanager dist
	
