CC = cc
CFLAGS ?= -O3 -std=gnu99
CPPFLAGS += -I. -DCLPS_NO_REJ_PRINTF
CORE = clps.c clps_packing.c fips202.c ntt.c packing.c poly.c polyvec.c reduce.c rounding.c symmetric-shake.c
HEADERS = $(wildcard *.h)
MODES = 2 3 5

.PHONY: all check clean
all: $(foreach m,$(MODES),build/bench-m$(m) build/bench-raw-m$(m) build/smoke-m$(m))

build:
	mkdir -p build

build/bench-m%: bench.c $(CORE) $(HEADERS) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDILITHIUM_MODE=$* bench.c $(CORE) -o $@

build/bench-raw-m%: bench_raw.c $(CORE) $(HEADERS) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDILITHIUM_MODE=$* bench_raw.c $(CORE) -o $@

build/smoke-m%: smoke.c $(CORE) $(HEADERS) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDILITHIUM_MODE=$* smoke.c $(CORE) -o $@

check: $(foreach m,$(MODES),build/smoke-m$(m))
	./build/smoke-m2
	./build/smoke-m3
	./build/smoke-m5

clean:
	rm -f build/bench-m2 build/bench-m3 build/bench-m5 build/bench-raw-m2 build/bench-raw-m3 build/bench-raw-m5 build/smoke-m2 build/smoke-m3 build/smoke-m5
