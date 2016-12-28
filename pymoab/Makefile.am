
EXTRA_DIST = $(srcdir)/cmake/cython_test.pyx $(srcdir)/pymoab/* $(srcdir)/tests/*


all-local:
	(cp --no-preserve=all -r $(srcdir)/pymoab .)
	(cd pymoab ; cython --cplus *.pyx)
	CFLAGS='-w -fPIC' python setup.py build
	(export PYTHONPATH=$(PWD):$(PYTHONPATH) ; python setup.py install --install-platlib=$(PWD))

install:
	if [ ! -d $(pymoab_prefix) ]; then mkdir -p $(pymoab_prefix); fi
	(export PYTHONPATH=$(pymoab_prefix):$(PYTHONPATH) ; python setup.py install --install-platlib=$(pymoab_prefix) )

check:
	(cp --no-preserve=all -r $(srcdir)/tests .)
	(cd tests ; PYTHONPATH=$(PWD):$(PYTHONPATH) nosetests)

uninstall:
	rm -r ${pymoab_prefix}
