
EXTRA_DIST = $(srcdir)/cmake/cython_test.pyx $(srcdir)/pymoab/* $(srcdir)/tests/*

local_install_path = local_inst

#CLEANFILES = pymoab/* pymoab build/* build tests/* tests $(local_install_path) pymoab.egg-info dist/* dist

test_files = test_core.py test_hom.py test_scd.py test_range.py test_tag.py test_topo_util.py

all-local:
	(if [ ! -d pymoab ]; then cp --no-preserve=all -r $(srcdir)/pymoab .; fi )
	(cd pymoab ; cython --cplus *.pyx)
	CFLAGS='-w -fPIC' python setup.py build
	(export PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) ; mkdir -m 766 $(local_install_path); python setup.py install --install-platlib=$(PWD)/$(local_install_path)/)

install:
	if [ ! -d $(pythondir) ]; then mkdir -m 766 -p $(pythondir); fi
	(export PYTHONPATH=$(pythondir):$(PYTHONPATH) ; python setup.py install --install-platlib=$(pythondir) )

check:
	(if [ ! -d tests ]; then cp --no-preserve=all -r $(srcdir)/tests .; fi)
	cd tests; \
	for i in $(test_files); \
	do \
		(LD_LIBRARY_PATH=$(abs_top_builddir)/src/.libs/:$(LD_LIBRARY_PATH) PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) python $$i) || exit 1; \
	done

clean-local:
	-rm -rf pymoab build tests $(local_install_path) pymoab.egg-info dist


uninstall:
	-rm -r ${pythondir}
