# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    	=
SPHINXBUILD   	= sphinx-build
SPHINXAPIDOC	= sphinx-apidoc
SOURCEDIR     	= source
MODULEDIR     	= ../../PythonPackage/
BUILDDIR      	= build
APIDOCOPTS  	= --force -M --separate
DOXYGENDIR		= ../doxygen
DOXYFILE		= DoxygenConfig
DOXYGEN			= doxygen

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	-cd $(DOXYGENDIR) && rm -rf build/xml build/html && doxygen DoxygenConfig && cp images/* build/xml/ || cp -r build/html ../sphinx/source/_static/cpp_doxygen_html
	rm -rf $(BUILDDIR)
	$(SPHINXAPIDOC) "$(MODULEDIR)" $(APIDOCOPTS) -o "$(SOURCEDIR)/_implementation" 
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
