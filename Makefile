# LMTsolver build
#
# Usage:
#   make                              # default: WernerAu, LMAX=51
#   make MATERIAL=DrudeAl
#   make MATERIAL=BrendelBormannAu LMAX=51
#   make LMAX=61
#   make all_materials                # build every material at current LMAX
#   make clean                        # remove all build artifacts
#   make info                         # print current settings
#
# Available MATERIAL values:
#   BrendelBormannAu  DrudeAl  JohnsonChristyAu  LorentzAu
#   LorentzAuv2  WernerAu  WernerAg  WernerBi
#
# Available LMAX values (must match a compiled integral table):
#   11  16  31  51  61  71

MATERIAL ?= WernerAu
LMAX     ?= 51

VALID_MATERIALS := DrudeAl LorentzSi BrendelBormannSi \
                   WernerAu WernerAg WernerBi \
                   BrendelBormannSiC OroscoCoimbraSiC 
VALID_LMAX      := 11 16 31 51 61 71

ifeq ($(filter $(MATERIAL),$(VALID_MATERIALS)),)
  $(error Unknown MATERIAL='$(MATERIAL)'. Valid: $(VALID_MATERIALS))
endif
ifeq ($(filter $(LMAX),$(VALID_LMAX)),)
  $(error Unknown LMAX='$(LMAX)'. Valid: $(VALID_LMAX))
endif

CXX      := g++
CXXFLAGS := -O2 -w -Iinclude -Iinclude/Dielectric/$(MATERIAL)
LDFLAGS  := -lcomplex_bessel

BUILDDIR := build
BINDIR   := bin
TARGET   := $(BINDIR)/$(MATERIAL)_LMTsolver

IN_OBJ   := $(BUILDDIR)/IN$(LMAX).o
MAIN_OBJ := $(BUILDDIR)/main_$(MATERIAL).o

# Materials that need Faddeeva
FADEVA_MATERIALS := BrendelBormannAu BrendelBormannSiC OroscoCoimbraSiC OroscoCoimbraAu

ifeq ($(filter $(MATERIAL),$(FADEVA_MATERIALS)),)
  EXTRA_OBJ :=
else
  EXTRA_OBJ := $(BUILDDIR)/Faddeeva.o
endif

OBJS := $(IN_OBJ) $(MAIN_OBJ) $(EXTRA_OBJ)

# ── Targets ─────────────────────────────────────────────────────────────────

.PHONY: all all_materials clean info

all: $(BUILDDIR) $(BINDIR) $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)
	@echo "Built: $@"

# main.cpp is compiled per-material (Epsilon.h is resolved via -I flag)
$(MAIN_OBJ): src/main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Integral table objects are material-independent
$(IN_OBJ): src/IN$(LMAX).cpp include/IN$(LMAX).h
	$(CXX) -Iinclude -w -O2 -c $< -o $@

$(BUILDDIR)/Faddeeva.o: src/Faddeeva.cpp include/Faddeeva.hh
	$(CXX) -Iinclude -w -O2 -c $< -o $@

$(BUILDDIR) $(BINDIR):
	mkdir -p $@

# Build every material at the chosen LMAX
all_materials:
	@for mat in $(VALID_MATERIALS); do \
	  echo "=== Building $$mat (LMAX=$(LMAX)) ==="; \
	  $(MAKE) --no-print-directory MATERIAL=$$mat LMAX=$(LMAX); \
	done

clean:
	rm -f $(BUILDDIR)/*.o $(BINDIR)/*

info:
	@echo "MATERIAL  : $(MATERIAL)"
	@echo "LMAX      : $(LMAX)"
	@echo "TARGET    : $(TARGET)"
	@echo "OBJECTS   : $(OBJS)"
