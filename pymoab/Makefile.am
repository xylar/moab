


all-local:
	(cp -r $(srcdir)/pymoab .)
	(cd pymoab ; cython --cplus *.pyx)
	CFLAGS='-w' python setup.py build
	(export PYTHONPATH=$(PWD):$(PYTHONPATH) ; python setup.py install --install-platlib=$(PWD))

install:
	if [ ! -d @prefix@/pymoab ]; then mkdir @prefix@/pymoab; fi
	(export PYTHONPATH=@prefix@/pymoab:$(PYTHONPATH) ; python setup.py install --install-platlib=@prefix@/pymoab/ )

check:
	(cp -r $(srcdir)/tests .)
	(PYTHONPATH=$(PWD):$(PYTHONPATH) cd tests ; nosetests)
