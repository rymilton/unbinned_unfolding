#===============================================================================
# Description:
#   Makefile for the RooUnfold package
#
# Instructions:
# o Make sure that ROOT (root-config actually) is installed in your PATH, or else
#   that the ROOTSYS environment variable is set and points to your ROOT release.
#   If not set up by default, this can be done
#   with:      . PATH_TO_ROOT/bin/thisroot.sh      (Bourne-style shell)
#   or with:   source PATH_TO_ROOT/bin/thisroot.csh     (C-style shell)
#
# o run "make TARGET"
#   - Default TARGET makes shared library (libRooUnfold.so) for use in ROOT.
#   - Use "make bin" to also build example programs.
#   - Use "make FILE" to compile and link user code in FILE.cxx .
#   - Add NOROOFIT=1 to link without RooFit (this is default if RooFit is not available).
#   - Add SHARED=1 to link test executables with shared library (libRooUnfold.so).
#     Otherwise links with static library (libRooUnfold.a).
#     To used the shared library: export LD_LIBRARY_PATH="$PWD:$LD_LIBRARY_PATH"
#   - Add ROOTBUILD=debug for debug version.
#   - Add VERBOSE=1 to show commands as they are executed.
#
# Build targets:
#   help    - give brief help
#   shlib   - make libRooUnfold.so (default target)
#   depend  - make dependency files (*.d)
#   lib     - make libRooUnfold.a
#   bin     - make lib and example programs
#   commands- show commands to make each type of target
#   html    - make documentation in htmldoc subdirectory
#   python  - make python package
#   cleanbin- delete test binaries and objects
#   clean   - delete all intermediate and final build objects
#   FILE    - if FILE.cxx (or FILE.cc or FILE.C) exists, builds executable FILE
#
# Author: Tim Adye <T.J.Adye@rl.ac.uk>
#
#===============================================================================

# === ROOT setup ===============================================================

RC           := root-config
ROOTSYS0     := $(shell $(RC) --prefix 2>/dev/null)
ifeq ($(ROOTSYS0),)
ifeq ($(ROOTSYS),)
$(error $(RC) not found and $$ROOTSYS is not defined)
endif
ETC          := $(ROOTSYS)/etc $(ROOTSYS)/test
else
ETC          := $(shell $(RC) --etcdir 2>/dev/null)
ifeq ($(ETC),)
ETC          := $(ROOTSYS0)/etc $(ROOTSYS0)/test
else
ETC          := $(ETC) $(ROOTSYS0)/test
endif
endif

ROOT_MAKEFILE=$(firstword $(wildcard $(addsuffix /Makefile.arch,$(ETC))))
ifeq ($(ROOT_MAKEFILE),)
ifeq ($(MAKE_RESTARTS),)
$(warning Could not find Makefile.arch under $(word 1,$(ETC)) or $(word 2,$(ETC)) - trying a basic Linux config)
endif
else
RC           :=
-include $(ROOT_MAKEFILE)
ifeq ($(ARCH),)
$(warning $(ROOT_MAKEFILE) setup failed (no $$ARCH) - trying a basic Linux config)
endif
endif
ifeq ($(RC),)
RC           := $(ROOTSYS)/bin/root-config
ifeq ($(wildcard $(RC)),)
$(error $(RC) not found)
endif
endif

ifeq ($(ARCH),)
# === This section is just in case ROOT's Makefile.arch didn't work =======
ARCH          = $(shell $(RC) --arch)
ROOTLIBS      = $(shell $(RC) --libs)
CXX           = $(shell $(RC) --cxx)
CXXFLAGS      = $(shell $(RC) --cflags)
CXXFLAGS     += -Wall -fPIC
LD            = $(shell $(RC) --ld)
LDFLAGS       = $(shell $(RC) --ldflags)
SOFLAGS       = -shared
ObjSuf        = o
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"
ifneq ($(findstring debug,$(ROOTBUILD)),)
CXXFLAGS     += -g
LDFLAGS      += -g
else
CXXFLAGS     += -O2
endif
endif

# add more warnings to the C++ compiler output
ifneq ($(findstring clang++,$(CXX)),)
# NB. don't include clang's -Wshadow, which gives many errors inside ROOT
CXXFLAGS     += -Wall -Wextra -Woverloaded-virtual
else
# NB. don't include gcc's -Wsuggest-override, which gives many warnings inside ROOT.
# Clang (and Cling)'s default and supposedly equivalent -Winconsistent-missing-override is more helpful.
CXXFLAGS     += -Wall -Wextra -Woverloaded-virtual -Wshadow
endif

ifeq ($(PLATFORM),macosx)
# Remove stupid shared library option on MacOSX. The option doesn't work (and we don't
# need it) because our linker output file specifies the full path.
SOFLAGS      := $(subst -install_name $(CURDIR)/,,$(SOFLAGS))
endif

ROOTINCDIR    = $(shell $(RC) --incdir)
ROOTLIBDIR    = $(shell $(RC) --libdir)
ROOTINCLUDES  = -I$(ROOTINCDIR)
ifeq ($(shell $(RC) --has-cling),yes)
ifeq ($(ROOTCLING),)
ROOTCLING     = rootcling
endif
else
ifeq ($(ROOTCINT),)
ROOTCINT      = rootcint
endif
ifeq ($(RLIBMAP),)
RLIBMAP       = rlibmap
endif
endif
ifeq ($(VERBOSE),1)
_             =
else
_             = @
endif

ifeq ($(MAKE_RESTARTS),)
$(info Use ROOT $(shell $(RC) --version) for $(ARCH) from $(shell $(RC) --prefix))
endif

# === RooUnfold directories ========================================

PACKAGE       = RooUnfold
SRCDIR        = $(CURDIR)/src/
INCDIR        = $(SRCDIR)
WORKDIR       = $(CURDIR)/build/$(ARCH)/
LIBDIR        = $(CURDIR)/
SHLIBDIR      = $(CURDIR)/
EXEDIR        =
EXESRC1       = $(CURDIR)/examples/
EXESRC2       = $(CURDIR)/test/src/
EXESRC        = $(EXESRC1)RooUnfoldExample.cxx $(EXESRC2)RooUnfoldTest.cxx $(EXESRC2)RooUnfoldTest2D.cxx $(EXESRC2)RooUnfoldTest3D.cxx
HTMLDOC       = htmldoc
OBJDIR        = $(WORKDIR)obj/
DEPDIR        = $(WORKDIR)dep/

# === RooUnfold options ===================================================

# Only provide RooUnfoldTUnfold if TUnfold is available.
# TUnfold is included in ROOT 5.22 and later.
ifeq ($(HAVE_TUNFOLD),)
ifneq ($(wildcard $(ROOTINCDIR)/TUnfold.h),)
HAVE_TUNFOLD  = 1
endif
endif

# Only use TUnfoldSys if it is available (ROOT 5.23.02 and later).
ifeq ($(HAVE_TUNFOLDSYS),)
ifneq ($(wildcard $(ROOTINCDIR)/TUnfoldSys.h),)
HAVE_TUNFOLDSYS = 1
endif
endif

ifneq ($(HAVE_TUNFOLD),1)
CPPFLAGS     += -DNOTUNFOLD
EXCLUDE      += RooUnfoldTUnfold.cxx RooUnfoldTUnfold.h
else
ifneq ($(HAVE_TUNFOLDSYS),1)
CPPFLAGS     += -DNOTUNFOLDSYS
endif
# Add -lUnfold if the library is available (added in ROOT 6.10)
ROOTLIBS     += $(patsubst $(ROOTLIBDIR)/lib%.$(DllSuf),-l%,$(wildcard $(patsubst %,$(ROOTLIBDIR)/lib%.$(DllSuf),Unfold)))
endif

# RooFit is included in ROOT if ROOT was compiled with --enable-roofit.
# We only use it for better-normalised test distributions in RooUnfoldTest
# (uses examples/RooUnfoldTestPdfRooFit.cxx instead of examples/RooUnfoldTestPdf.cxx).
ifeq ($(NOROOFIT),)
ifneq ($(shell $(RC) --has-roofit),yes)
$(warning This version of ROOT does not support RooFit. We will build the test programs without it.)
NOROOFIT      = 1
endif
endif

ifneq ($(NOROOFIT),)
CPPFLAGS     += -DNOROOFIT
else
# Different versions of ROOT require different libraries with RooFit, so
# include all the ones that exist.
# Note that if the ROOT shared libraries were linked against them
# (configure --enable-explicitlink ?), as is done in the CERN AFS versions,
# then these are not required. But they do no harm.
ROOFITLIBS   += -lRooFit -lRooFitCore -lHistFactory $(patsubst $(ROOTLIBDIR)/lib%.$(DllSuf),-l%,$(wildcard $(patsubst %,$(ROOTLIBDIR)/lib%.$(DllSuf),Thread Minuit Foam MathMore Html)))
endif

# === Internal configuration ===================================================

MAIN          = $(notdir $(EXESRC))
MAINEXE       = $(addprefix $(EXEDIR),$(patsubst %.cxx,%$(ExeSuf),$(MAIN)))
LINKDEF       = $(INCDIR)$(PACKAGE)_LinkDef.h
LINKDEFMAP    = $(WORKDIR)$(PACKAGE)Map_LinkDef
HLIST         = $(filter-out $(addprefix $(INCDIR),$(EXCLUDE)) $(LINKDEF),$(wildcard $(INCDIR)*.h $(INCDIR)*.tpp)) $(LINKDEF)
CINTFILE      = $(WORKDIR)$(PACKAGE)Dict.cxx
ifneq ($(ROOTCLING),)
CLINGDICT     = $(SHLIBDIR)$(PACKAGE)Dict_rdict.pcm
endif
CINTOBJ       = $(OBJDIR)$(PACKAGE)Dict.$(ObjSuf)
LIBNAME       = $(PACKAGE)
STATICLIBNAME = $(PACKAGE)_static
LIBFILE       = $(LIBDIR)lib$(LIBNAME).a
STATICLIBFILE = $(LIBDIR)lib$(STATICLIBNAME).a
SHLIBFILE     = $(SHLIBDIR)lib$(LIBNAME).$(DllSuf)
ROOTMAP       = $(SHLIBDIR)lib$(LIBNAME).rootmap
ifeq ($(MAKECMDGOALS),)
GOALS         = default
else
GOALS         = $(filter-out clean cleanbin help commands,$(MAKECMDGOALS))
endif

CPPFLAGS     += -DMAKEBUILD
ifneq ($(findstring g++,$(CXX)),)
MFLAGS        = -MM
endif
INCLUDES      = -I$(INCDIR) -I$(CURDIR)
CXX          += $(EXTRAINCLUDES)
LDFLAGS      += $(EXTRALDFLAGS)

ifneq ($(SHARED),)
LIBS          = -L$(SHLIBDIR)
LINKLIB       = $(SHLIBFILE)
LINKLIBOPT    = -l$(LIBNAME) $(EXTRALIBS)
else
LIBS          = -L$(LIBDIR)
LINKLIB       = $(STATICLIBFILE)
LINKLIBOPT    = -l$(STATICLIBNAME) $(EXTRALIBS)
endif

# List of all object files to build
SRCLIST       = $(filter-out $(EXCLUDE),$(notdir $(wildcard $(SRCDIR)*.cxx))) $(EXTRASRC)
OLIST         = $(addprefix $(OBJDIR),$(addsuffix .$(ObjSuf),$(basename $(SRCLIST))))

ifneq ($(filter %.for,$(SRCLIST)),)
  ifeq ($(findstring 77,$(FC)),77)
	GCCLIBS       = -lg2c
  endif	
endif

ROOTLIBFILES := $(wildcard $(patsubst -l%,$(ROOTLIBDIR)/lib%.$(DllSuf),$(filter -l%,$(ROOTLIBS))))

ifeq ($(MFLAGS),)

# Can't make dependency files, so make every compilation dependent on all headers.
HDEP          = $(HLIST)

else

# List of all dependency files to make
DLIST         = $(addprefix $(DEPDIR),$(patsubst %.cxx,%.d,$(filter %.cxx,$(SRCLIST) $(MAIN))))

# If possible, limit ROOTLIBFILES to libraries that we actually use.
DLISTLIB      = $(wildcard $(addprefix $(DEPDIR),$(patsubst %.cxx,%.d,$(filter %.cxx,$(SRCLIST)))))
ifneq ($(DLISTLIB),)
ROOTLIBMAPS  := $(wildcard $(patsubst %.$(DllSuf),%.rootmap,$(ROOTLIBFILES)))
ifneq ($(ROOTLIBMAPS),)
ROOTLIBFILESUSED := $(patsubst %.rootmap,%.$(DllSuf),$(shell grep -l '^Library\.\($(shell sed -n 's,^ *$(ROOTINCDIR)[^ ]*/\([^ /]*\)\.h .*,\1\\,p' $(DLISTLIB) | sort -u | tr '\n' '|' | sed 's,\\|$$,,')\):' $(ROOTLIBMAPS)))
ifneq ($(ROOTLIBFILESUSED),)
ROOTLIBFILES := $(ROOTLIBFILESUSED)
endif
endif
endif

endif

ifeq ($(NOROOFIT),)
ifeq ($(wildcard $(DLIST)),)
# Since we can't check the dependencies, include RooFit on all binaries
ROOFITCLIENTS = $(patsubst %.cxx,$(OBJDIR)%.$(ObjSuf),$(MAIN))
else
# List of programs that use RooFit. Should only be those in $(EXESRC)
ROOFITCLIENTS = $(patsubst $(DEPDIR)%.d,$(OBJDIR)%.$(ObjSuf),$(shell fgrep -l '/RooAbsArg.h ' $(DLIST) 2>/dev/null))
endif
endif

# === Implicit rules ===========================================================

# Remove default compile-and-link rules
%$(ExeSuf)  : %.cc
%$(ExeSuf)  : %.C
%$(ExeSuf)  : %.$(SrcSuf)
%$(ExeSuf)  : %.$(ObjSuf)

# Implicit rule making all dependency Makefiles included at the end of this makefile
$(DEPDIR)%.d : $(SRCDIR)%.cxx
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@rm -f $@
	$(_)set -e; \
	 $(CXX) $(MFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< \
	 | sed 's,\($(notdir $*)\.$(ObjSuf)\) *:,$(OBJDIR)\1 $@ :,g' > $@; \
	 [ -s $@ ] || rm -f $@

$(DEPDIR)%.d : $(EXESRC1)%.cxx
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@rm -f $@
	$(_)set -e; \
	 $(CXX) $(MFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< \
	 | sed 's,\($(notdir $*)\.$(ObjSuf)\) *:,$(OBJDIR)\1 $@ :,g' > $@; \
	 [ -s $@ ] || rm -f $@

$(DEPDIR)%.d : $(EXESRC2)%.cxx
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@rm -f $@
	$(_)set -e; \
	 $(CXX) $(MFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< \
	 | sed 's,\($(notdir $*)\.$(ObjSuf)\) *:,$(OBJDIR)\1 $@ :,g' > $@; \
	 [ -s $@ ] || rm -f $@

# Implicit rule to compile all classes
$(OBJDIR)%.$(ObjSuf) : $(SRCDIR)%.cxx $(HDEP)
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(INCLUDES)

# Implicit rule to compile intermediate files
$(OBJDIR)%.$(ObjSuf) : $(WORKDIR)%.cxx $(HDEP)
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(INCLUDES)

# Implicit rule for Fortran (if any)
$(OBJDIR)%.$(ObjSuf) : $(SRCDIR)%.for $(FDEP)
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	$(_)$(FC) $(FFFLAGS) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Implicit rule to compile main program
$(OBJDIR)%.$(ObjSuf) : $(EXESRC1)%.cxx $(HLIST)
	@echo "Compiling example program $<"
	@mkdir -p $(OBJDIR)
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(INCLUDES)

$(OBJDIR)%.$(ObjSuf) : $(EXESRC2)%.cxx $(HLIST)
	@echo "Compiling example program $<"
	@mkdir -p $(OBJDIR)
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(INCLUDES)

%.$(ObjSuf) : %.$(SrcSuf) $(HLIST)
	@echo "Compiling user program $<"
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(INCLUDES)

%.$(ObjSuf) : %.cc $(HLIST)
	@echo "Compiling user program $<"
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(INCLUDES)

%.$(ObjSuf) : %.C $(HLIST)
	@echo "Compiling user program $<"
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(INCLUDES)

# Implicit rules to link main program
$(EXEDIR)%$(ExeSuf) : $(OBJDIR)%.$(ObjSuf) $(LINKLIB)
	@echo "Linking example executable $@"
	$(_)$(LD) $(LDFLAGS) $< $(OutPutOpt)$@ $(LIBS) $(LINKLIBOPT) $(ROOTLIBS) $(if $(findstring $<,$(ROOFITCLIENTS)),$(ROOFITLIBS)) $(GCCLIBS)

%$(ExeSuf) : %.$(ObjSuf) $(LINKLIB)
	@echo "Linking user executable $@"
	$(_)$(LD) $(LDFLAGS) $< $(OutPutOpt)$@ $(LIBS) $(LINKLIBOPT) $(ROOTLIBS) $(if $(findstring $<,$(ROOFITCLIENTS)),$(ROOFITLIBS)) $(GCCLIBS)

# === Explicit rules ===========================================================

default : shlib python

help        :
	@echo "Usage: $(MAKE) [TARGET] [ROOTBUILD=debug] [VERBOSE=1] [NOROOFIT=1] [SHARED=1]"
	@echo "Some TARGETs are: 'bin', 'html', 'clean', 'python', and 'commands'"

# Rule to make ROOTCINT output file
ifeq ($(ROOTCLING),)
SHLIBEXTRA=$(ROOTMAP)
$(CINTFILE) : $(HLIST)
	@mkdir -p $(WORKDIR)
	@mkdir -p $(OBJDIR)
	@echo "Generating CINT dictionary $(CINTFILE)"
	$(_)$(ROOTCINT) -f $(CINTFILE) -c -p $(CPPFLAGS) $(INCLUDES) $(HLIST)

$(ROOTMAP) : $(SHLIBFILE) $(LINKDEF)
	@echo "Making CINT $@"
	$(_)echo TObject.h TMemberInspector.h $(HLIST) | tr ' ' '\n' | sed 's/^\(.*\)$$/#include "\1"/' > $(LINKDEFMAP).cxx
	$(_)$(CXX) -E -D__MAKECINT__ -D__CINT__ -o $(LINKDEFMAP).h $(CPPFLAGS) $(INCLUDES) $(ROOTINCLUDES) $(LINKDEFMAP).cxx
	$(_)$(RLIBMAP) -o $@ -l $< -d $(ROOTLIBFILES) -c $(LINKDEFMAP).h
else
SHLIBEXTRA=
$(CINTFILE) : $(HLIST)
	@mkdir -p $(WORKDIR)
	@mkdir -p $(OBJDIR)
	@echo "Generating CLING dictionary files $(CINTFILE), $(CLINGDICT), and $(ROOTMAP)"
	$(_)$(ROOTCLING) -f $(CINTFILE) -rml $(notdir $(SHLIBFILE)) -rmf $(ROOTMAP) $(CPPFLAGS) $(INCLUDES) $(HLIST)
	$(_)mv -f $(dir $(CINTFILE))$(notdir $(CLINGDICT)) $(CLINGDICT)
endif

# Rule to combine objects into a library
$(LIBFILE) : $(OLIST) $(CINTOBJ)
	@echo "Making $@"
	@mkdir -p $(LIBDIR)
	@rm -f $@
	$(_)ar cq $@ $^
	$(_)ranlib $@

# Make symlink to static library so we can refer to that without picking up shared library.
# Assumes $(LIBFILE) and $(STATICLIBFILE) are in the same directory
$(STATICLIBFILE) : $(LIBFILE)
	@rm -f $@
	$(_)ln -s $(patsubst $(dir $@)%,%,$<) $@

# Rule to combine objects into a shared library
$(SHLIBFILE) : $(OLIST) $(CINTOBJ)
	@echo "Making $@"
	@mkdir -p $(SHLIBDIR)
	@rm -f $@
	$(_)$(LD) $(SOFLAGS) $(ROOFITLIBS) $(LDFLAGS) $^ $(OutPutOpt)$@ $(ROOTLIBS) $(GCCLIBS)

# Useful build targets
depend: $(DLIST)
include: depend
lib: $(LIBFILE)
shlib: $(SHLIBFILE) $(SHLIBEXTRA)
exe: $(MAINEXE)
bin: shlib exe

commands :
	@echo "Make $(DEPDIR)%.d:	$(CXX) $(MFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $(SRCDIR)%.cxx | sed 's,\(%\.$(ObjSuf)\) *:,$(OBJDIR)\1 $(DEPDIR)%.d :,g' > $(DEPDIR)%.d"
	@echo
ifeq ($(ROOTCLING),)
	@echo "Make dictionary:         $(ROOTCINT) -f $(CINTFILE) -c -p $(CPPFLAGS) $(INCLUDES) $(HLIST)"
else
	@echo "Make dictionary:         $(ROOTCLING) -f $(CINTFILE) -rml $(notdir $(SHLIBFILE)) -rmf $(ROOTMAP) $(CPPFLAGS) $(INCLUDES) $(HLIST)"
endif
	@echo
	@echo "Compile $(SRCDIR)%.cxx:	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(SRCDIR)%.cxx -o $(OBJDIR)%.$(ObjSuf) $(INCLUDES)"
	@echo
	@echo "Make $(SHLIBFILE):	$(LD) $(SOFLAGS) $(LDFLAGS) *.$(ObjSuf) $(OutPutOpt)$(SHLIBFILE) $(ROOTLIBS) $(GCCLIBS)"
	@echo
	@echo "Make executable $(EXEDIR)RooUnfoldTest$(ExeSuf):	$(LD) $(LDFLAGS) $(OBJDIR)RooUnfoldTest.$(ObjSuf) $(OutPutOpt)$(EXEDIR)RooUnfoldTest$(ExeSuf) $(LIBS) $(LINKLIBOPT) $(ROOTLIBS) $(ROOFITLIBS) $(GCCLIBS)"

clean : cleanbin
	rm -f $(DLIST)
	rm -f $(CINTFILE) $(basename $(CINTFILE)).h $(CLINGDICT)
	rm -f $(OLIST) $(CINTOBJ)
	rm -f $(LIBFILE)
	rm -f $(SHLIBFILE) $(ROOTMAP) $(LINKDEFMAP).cxx $(LINKDEFMAP).h
	rm -f $(STATICLIBFILE)
	rm -f $(LIBDIR)$(PACKAGE)/__init__.py
	rm -fd $(LIBDIR)$(PACKAGE) $(OBJDIR) $(DEPDIR) $(WORKDIR)

cleanbin :
	rm -f $(addprefix $(OBJDIR),$(patsubst %.cxx,%.$(ObjSuf),$(MAIN)))
	rm -f $(MAINEXE)

$(HTMLDOC)/index.html : $(SHLIBFILE)
	@echo "Making HTML documentation in $(HTMLDOC)"
	@( echo 'gSystem->Load("lib$(LIBNAME)");'; \
	   echo 'THtml h;';                        \
	   echo 'h.SetProductName("$(PACKAGE)");'; \
	   echo 'h.SetInputDir("$(SRCDIR)");';     \
	   echo 'h.SetOutputDir("$(HTMLDOC)");';   \
	   echo 'h.MakeAll();';                    \
	   echo '.q' )                             \
	 | ROOT_HIST=0 root -n -l -b

html : $(HTMLDOC)/index.html

python : $(LIBDIR)/$(PACKAGE)/__init__.py

$(LIBDIR)/$(PACKAGE)/%.py: $(LIBDIR)/python/%.py
	@mkdir -p $(LIBDIR)/$(PACKAGE)
	@ln -nsf $^ $@

.PHONY : include depend shlib lib exe bin default clean cleanbin html help python

ifneq ($(GOALS),)
ifneq ($(DLIST),)
-include $(DLIST)
endif
endif
