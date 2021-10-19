all :
	python setup.py install --user

html:
	make -C doc html &> /dev/null
