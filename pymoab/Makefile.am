


all-local:
	(cd $(srcdir)/pymoab ; cython --cplus *.pyx)
	(cd $(srcdir); cp setup.py.in setup.py)
	(cd $(srcdir); sed -i 's,$${CMAKE_SOURCE_DIR},$(abspath $(top_srcdir)),g' setup.py)
	(cd $(srcdir); sed -i 's,$${CMAKE_BINARY_DIR}/lib/,$(abspath $(top_builddir))/src/.libs/,g' setup.py)
	(cd $(srcdir); sed -i 's,$${CMAKE_BINARY_DIR},$(abspath $(top_builddir)),g' setup.py)
	(cd $(srcdir); sed -i 's,$${CMAKE_CURRENT_BINARY_DIR},$(abspath $(srcdir)),g' setup.py)
	(cd $(srcdir); python setup.py build \
	--build-base $(shell readlink -f $(builddir))/build \
	--verbose)

