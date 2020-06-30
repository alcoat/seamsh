all :
	python setup.py install --user
	make -C doc html

html:
	make -C doc html
