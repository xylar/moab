


all-local:
	(cp -r $(srcdir)/pymoab .)
	(cd pymoab ; cython --cplus *.pyx)
	CFLAGS='-w' python setup.py build
	(export PYTHONPATH=$(PWD):$(PYTHONPATH) ; python setup.py install --install-platlib=$(PWD))

install:
	if [ ! -d $(pymoab_prefix) ]; then mkdir -p $(pymoab_prefix); fi
	(export PYTHONPATH=$(pymoab_prefix):$(PYTHONPATH) ; python setup.py install --install-platlib=$(pymoab_prefix) )

check:
	(cp -r $(srcdir)/tests .)
	(cd tests ; PYTHONPATH=$(PWD):$(PYTHONPATH) nosetests)
