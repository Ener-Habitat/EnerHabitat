# Standalone build of the 2D free-running C solver (no PostgreSQL, no mail).
# Golden master / time baseline for the Python port (see PLAN-2D.md, Fase 0).
#
# readline/ncurses are pulled from Homebrew on macOS; override READLINE_PREFIX
# if installed elsewhere (e.g. /usr/local on Intel macs, /usr on Linux).

READLINE_PREFIX ?= $(shell brew --prefix readline 2>/dev/null)
CC      ?= gcc
# -std=gnu89: the legacy code uses K&R implicit-int return types, which modern
# clang rejects by default. gnu89 keeps them as warnings.
CFLAGS  ?= -O2 -std=gnu89 -Wno-implicit-int -Wno-implicit-function-declaration
INCS     = $(if $(READLINE_PREFIX),-I$(READLINE_PREFIX)/include)
LIBS     = $(if $(READLINE_PREFIX),-L$(READLINE_PREFIX)/lib) -lreadline -lncurses -lm
DEFS     = -DREADLINE -DREUSE -DSOURCE=\"$<\" -DDATE=\""dev"\" -DHOSTNAME=\"standalone\"

# Build with `make` for the plain solver, or `make DUMP=1` to also write the
# geometry/material dumps (dat/dump_*.dat) used to verify the Python port.
ifdef DUMP
DEFS += -DDUMP
endif

# `make DUMPCOEF=1` dumps the coefficient assembly (a,b,c,d) for a known T/To
# field and exits immediately (no time loop) — for Fase 2 verification.
ifdef DUMPCOEF
DEFS += -DDUMPCOEF
endif

# `make DUMPSTEP=1` runs ONE solve_PQ step from a known field and exits — Fase 3.
ifdef DUMPSTEP
DEFS += -DDUMPSTEP
endif

# Day-to-day convergence cap. Default: converge to periodic steady state.
# `make DAYS=1` runs a single day; `make DAYS=N` runs N days.
ifdef DAYS
DEFS += -DMAXDAYS=$(DAYS)
endif

conduction.e: conduction.c arrays.h sol.h new_input.h tools.h | dat
	$(CC) $< -o $@ $(CFLAGS) $(DEFS) $(INCS) $(LIBS)

dat:
	mkdir -p dat

.PHONY: clean
clean:
	rm -f conduction.e
