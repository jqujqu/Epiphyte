#    Copyright (C) 2015 Andrew D. Smith
#
#    Authors: Jenny Qu, Andrew D. Smith
#
#    This code is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(realpath ../smithlab_cpp)
endif

ifndef TREETOOL
TREETOOL=$(realpath ../adssrc/treetool)
endif

PROGS = phylo-epi-tester simulation indep-sites-hypo

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50
OPTFLAGS = -O2 
DEBUGFLAGS = -g
LIBS = -lgsl -lgslcblas

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

COMMON_DIR = $(realpath ../common)
INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_DIR) ${TREETOOL}
INCLUDEARGS = $(addprefix -I, $(INCLUDEDIRS))

# Flags for Mavericks
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif

all: $(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o)  

simulation: $(addprefix $(COMMON_DIR)/, PhyloTreePreorder.o ) \
	$(addprefix $(TREETOOL)/, PhyloTree.o )

phylo-epi-tester: $(addprefix $(COMMON_DIR)/, MethpipeFiles.o BetaBin.o PhyloTreePreorder.o) \
	$(addprefix $(TREETOOL)/, PhyloTree.o )

indep-sites-hypo: $(addprefix $(COMMON_DIR)/, MethpipeFiles.o BetaBin.o) 

install: all
	@mkdir -p ../../bin
	@install -m 755 $(PROGS) ../../bin

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
