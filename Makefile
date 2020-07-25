all :
	make -C seamshlib/build
	cp seamshlib/build/libseamsh.so seamsh/
	python setup.py install --user

html:
	make -C doc html &> /dev/null
