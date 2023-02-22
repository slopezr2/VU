#
# Makefile for main directory.
#

help:
	@echo " "
	@echo " Usage:"
	@echo "   make docu           # build documentation"
	@echo "   make clean          # remove temporary files and generated documentation"
	@echo " "

# generate documentation:
docu:
	@echo " "
	(cd doc; make html)
	@echo " "
	@echo "Browse through documentation:"
	@echo "   file://$(PWD)/doc/build/html/index.html"
	@echo " "

# remove temporary files:
clean: clean-pyc clean-osx clean-docu
	/bin/rm -f setup_le*.*.log
	/bin/rm -f jb.[0-9]*.*
	/bin/rm -f lekf.copy.jb
	/bin/rm -f tmp-ttb-*.rc
	/bin/rm -f tmp-fc-*.rc

# remove temporary python objects
clean-pyc:
	rm -f -r `find . -name '__pycache__'`
	rm -f `find . -name '*.pyc'`

# remove temporary Mac files
clean-osx:
	rm -f `find . -name '.DS_Store'`

# remove temporary and build documentation files:
clean-docu:
	rm -f -r doc/build/doctrees
	rm -f -r doc/build/html

