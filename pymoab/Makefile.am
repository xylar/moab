AUTOMAKE_OPTIONS=foreign

EXTRA_DIST = $(srcdir)/cmake/cython_test.pyx $(srcdir)/pymoab/* $(srcdir)/tests/*

local_install_path = local_inst

#CLEANFILES = pymoab/* pymoab build/* build tests/* tests $(local_install_path) pymoab.egg-info dist/* dist

test_files = test_core.py test_hom.py test_scd.py test_range.py test_tag.py test_topo_util.py test_skinner.py

if DEVELOPER_MODE
all-local: cythonize
else
all-local:
endif
	CFLAGS='-w -fPIC' $(PYTHON) setup.py build
	($(MKDIR_P) -m 766 $(local_install_path); \
		PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) \
		$(PYTHON) setup.py install --install-platlib=$(PWD)/$(local_install_path)/ --install-purelib=$(PWD)/$(local_install_path) \
	)

install-exec-local: all-local
	$(MKDIR_P) $(DESTDIR)$(pkgpythondir)
	(PYTHONPATH=$(DESTDIR)$(pkgpythondir):$(PYTHONPATH) $(PYTHON) setup.py install \
		--prefix=$(DESTDIR)$(prefix) \
		--single-version-externally-managed \
		--record $(DESTDIR)$(pkgpythondir)/install_files.txt \
		--verbose \
	)

if DEVELOPER_MODE

$(srcdir)/pymoab/core.cpp: $(srcdir)/pymoab/core.pyx $(srcdir)/pymoab/moab.pxd
	$(CYTHON) --cplus --fast-fail --include-dir $(top_srcdir)/src --output-file $@ $<
$(srcdir)/pymoab/hcoord.cpp: $(srcdir)/pymoab/hcoord.pyx $(srcdir)/pymoab/moab.pxd
	$(CYTHON) --cplus --fast-fail --include-dir $(top_srcdir)/src --output-file $@ $<
$(srcdir)/pymoab/rng.cpp: $(srcdir)/pymoab/rng.pyx $(srcdir)/pymoab/moab.pxd
	$(CYTHON) --cplus --fast-fail --include-dir $(top_srcdir)/src --output-file $@ $<
$(srcdir)/pymoab/scd.cpp: $(srcdir)/pymoab/scd.pyx $(srcdir)/pymoab/moab.pxd
	$(CYTHON) --cplus --fast-fail --include-dir $(top_srcdir)/src --output-file $@ $<
$(srcdir)/pymoab/tag.cpp: $(srcdir)/pymoab/tag.pyx $(srcdir)/pymoab/moab.pxd
	$(CYTHON) --cplus --fast-fail --include-dir $(top_srcdir)/src --output-file $@ $<
$(srcdir)/pymoab/topo_util.cpp: $(srcdir)/pymoab/topo_util.pyx $(srcdir)/pymoab/moab.pxd
	$(CYTHON) --cplus --fast-fail --include-dir $(top_srcdir)/src --output-file $@ $<
$(srcdir)/pymoab/types.cpp: $(srcdir)/pymoab/types.pyx $(srcdir)/pymoab/moab.pxd
	$(CYTHON) --cplus --fast-fail --include-dir $(top_srcdir)/src --output-file $@ $<
$(srcdir)/pymoab/skinner.cpp: $(srcdir)/pymoab/skinner.pyx $(srcdir)/pymoab/moab.pxd
	$(CYTHON) --cplus --fast-fail --include-dir $(top_srcdir)/src --output-file $@ $<

cythonize: $(srcdir)/pymoab/core.cpp $(srcdir)/pymoab/hcoord.cpp \
					 $(srcdir)/pymoab/rng.cpp $(srcdir)/pymoab/scd.cpp \
					 $(srcdir)/pymoab/tag.cpp $(srcdir)/pymoab/topo_util.cpp \
					 $(srcdir)/pymoab/types.cpp $(srcdir)/pymoab/skinner.cpp

endif

check: all-local
	cd tests; \
	for i in $(test_files); \
	do \
		(@SLIBENVPATH@=$(abs_top_builddir)/src/.libs/:$(@SLIBENVPATH@) PYTHONPATH=$(abs_builddir)/$(local_install_path)/:$(PYTHONPATH) $(PYTHON) $$i) || exit 1; \
	done

clean-local:
	-rm -rf build tests/*.pyc tests/outfile.h5m $(local_install_path) pymoab.egg-info dist

uninstall:
	cat $(DESTDIR)$(pkgpythondir)/install_files.txt | xargs rm -rf
	-rm -rf $(DESTDIR)$(pkgpythondir)

