#
# Makefile for TM5 source directory.
#

help:
	@echo " "
	@echo " Usage:"
	@echo "   make clean          # remove temporary jb.* files"
	@echo " "

clean: clean-pyc
	/bin/rm -f setup_le*.*.log
	/bin/rm -f jb.[0-9]*.*
	/bin/rm -f tmp-ttb-*.rc
	/bin/rm -f tmp-fc-*.rc

clean-pyc:
	rm -f -r `find . -name '__pycache__'`
	rm -f `find . -name '*.pyc'`


