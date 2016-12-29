
EXTRA_DIST = $(srcdir)/cmake/cython_test.pyx $(srcdir)/pymoab/* $(srcdir)/tests/*

local_install_path = local_inst

all-local:
	(cp --no-preserve=all -r $(srcdir)/pymoab .)
	(cd pymoab ; cython --cplus *.pyx)
	CFLAGS='-w -fPIC' python setup.py build
	(export PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) ; mkdir -m 766 $(local_install_path); python setup.py install --install-platlib=$(PWD)/$(local_install_path)/)

install:
	if [ ! -d $(pythondir) ]; then mkdir -m 766 -p $(pythondir); fi
	(export PYTHONPATH=$(pythondir):$(PYTHONPATH) ; python setup.py install --install-platlib=$(pythondir) )

check:
	(cp --no-preserve=all -r $(srcdir)/tests .)
	(cd tests ; PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) nosetests)

clean-local:
	rm -r pymoab build tests $(local_install_path) pymoab.egg-info dist


uninstall:
	rm -r ${pythondir}
