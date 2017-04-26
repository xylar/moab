
EXTRA_DIST = $(srcdir)/cmake/cython_test.pyx $(srcdir)/pymoab/* $(srcdir)/tests/*

local_install_path = local_inst

#CLEANFILES = pymoab/* pymoab build/* build tests/* tests $(local_install_path) pymoab.egg-info dist/* dist

test_files = test_core.py test_hom.py test_scd.py test_range.py test_tag.py test_topo_util.py

all-local:
	CFLAGS='-w -fPIC' $(PYTHON) setup.py build
	(export PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) ; mkdir -m 766 $(local_install_path); $(PYTHON) setup.py install --install-platlib=$(PWD)/$(local_install_path)/ --install-purelib=$(PWD)/$(local_install_path)/)

install:
	if [ ! -d $(pythondir) ]; then mkdir -m 766 -p $(pythondir); fi
	(export PYTHONPATH=$(pythondir):$(PYTHONPATH) ; $(PYTHON) setup.py install --install-platlib=$(DESTDIR)$(pythondir) --install-purelib=$(DESTDIR)$(pythondir))

check:
	cd tests; \
	for i in $(test_files); \
	do \
		(LD_LIBRARY_PATH=$(abs_top_builddir)/src/.libs/:$(LD_LIBRARY_PATH) PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) $(PYTHON) $$i) || exit 1; \
	done

clean-local:
	-rm -rf build tests/*.pyc $(local_install_path) pymoab.egg-info dist

uninstall:
	-rm -r ${pythondir}
