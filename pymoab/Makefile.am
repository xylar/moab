
EXTRA_DIST = $(srcdir)/cmake/cython_test.pyx $(srcdir)/pymoab/* $(srcdir)/tests/*

local_install_path = local_inst

#CLEANFILES = pymoab/* pymoab build/* build tests/* tests $(local_install_path) pymoab.egg-info dist/* dist

test_files = test_core.py test_hom.py test_scd.py test_range.py test_tag.py test_topo_util.py test_skinner.py

all-local:
	CFLAGS='-w -fPIC' $(PYTHON) setup.py build
	(export PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) ; $(MKDIR_P) -m 766 $(local_install_path); $(PYTHON) setup.py install --install-platlib=$(PWD)/$(local_install_path)/ --install-purelib=$(PWD)/$(local_install_path)/)

install:
	$(MKDIR_P) $(DESTDIR)$(pkgpythondir)
	(export PYTHONPATH=$(DESTDIR)$(pkgpythondir):$(PYTHONPATH) ; $(PYTHON) setup.py install \
		--prefix=$(DESTDIR)$(prefix) \
		--single-version-externally-managed \
		--record $(DESTDIR)$(pkgpythondir)/install_files.txt \
		--verbose \
	)

%.cpp.pyx:
	cython --cplus -t -o $@ $<
 
$(srcdir)/pymoab/core.cpp: $(srcdir)/pymoab/core.pyx $(srcdir)/pymoab/moab.pxd
	cython --cplus -t -o $@ $<
$(srcdir)/pymoab/hcoord.cpp: $(srcdir)/pymoab/hcoord.pyx $(srcdir)/pymoab/moab.pxd
	cython --cplus -t -o $@ $<
$(srcdir)/pymoab/rng.cpp: $(srcdir)/pymoab/rng.pyx $(srcdir)/pymoab/moab.pxd
	cython --cplus -t -o $@ $<
$(srcdir)/pymoab/scd.cpp: $(srcdir)/pymoab/scd.pyx $(srcdir)/pymoab/moab.pxd
	cython --cplus -t -o $@ $<
$(srcdir)/pymoab/tag.cpp: $(srcdir)/pymoab/tag.pyx $(srcdir)/pymoab/moab.pxd
	cython --cplus -t -o $@ $<
$(srcdir)/pymoab/topo_util.cpp: $(srcdir)/pymoab/topo_util.pyx $(srcdir)/pymoab/moab.pxd
	cython --cplus -t -o $@ $<
$(srcdir)/pymoab/types.cpp: $(srcdir)/pymoab/types.pyx $(srcdir)/pymoab/moab.pxd
	cython --cplus -t -o $@ $<
$(srcdir)/pymoab/skinner.cpp: $(srcdir)/pymoab/skinner.pyx $(srcdir)/pymoab/moab.pxd
	cython --cplus -t -o $@ $<

cythonize: $(srcdir)/pymoab/core.cpp $(srcdir)/pymoab/hcoord.cpp $(srcdir)/pymoab/rng.cpp $(srcdir)/pymoab/scd.cpp $(srcdir)/pymoab/tag.cpp $(srcdir)/pymoab/topo_util.cpp $(srcdir)/pymoab/types.cpp $(srcdir)/pymoab/skinner.cpp

check:
	cd tests; \
	for i in $(test_files); \
	do \
		(LD_LIBRARY_PATH=$(abs_top_builddir)/src/.libs/:$(LD_LIBRARY_PATH) PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) $(PYTHON) $$i) || exit 1; \
	done

clean-local:
	-rm -rf build tests/*.pyc tests/outfile.h5m $(local_install_path) pymoab.egg-info dist

uninstall:
	cat $(DESTDIR)$(pkgpythondir)/install_files.txt | xargs rm -rf
	-rm -rf $(DESTDIR)$(pkgpythondir)

