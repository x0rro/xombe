BUILD ?= build
SAN_BUILD ?= build-asan
Q15_BUILD ?= build-q15
CMAKE ?= cmake
PYTHON ?= python3
CMAKE_FLAGS ?=
# WERROR=ON turns warnings into errors in every build
WERROR ?= OFF
JOBS ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

.PHONY: all check check-q15 conformance conformance-dstar clean

all:
	$(CMAKE) -S . -B $(BUILD) -DCMAKE_BUILD_TYPE=Release -DXOMBE_WERROR=$(WERROR) $(CMAKE_FLAGS)
	$(CMAKE) --build $(BUILD) -j$(JOBS)

# release, q15 (check-q15), and ASan/UBSan builds
check: all check-q15
	cd $(BUILD) && ctest --output-on-failure -j$(JOBS)
	$(CMAKE) -S . -B $(SAN_BUILD) -DCMAKE_BUILD_TYPE=Debug -DXOMBE_SANITIZE=ON -DXOMBE_WERROR=$(WERROR) $(CMAKE_FLAGS)
	$(CMAKE) --build $(SAN_BUILD) -j$(JOBS)
	cd $(SAN_BUILD) && ctest --output-on-failure -j$(JOBS)
	@n=$$(cd $(BUILD) && ctest -N -L no_python | sed -n 's/^Total Tests: //p'); \
	if [ "$${n:-0}" -gt 0 ]; then echo "make check: $$n tests skipped per build (no Python 3)"; fi

# the q15 kernels, with q15ops_arm when arm-none-eabi-gcc and qemu-arm are installed
check-q15:
	$(CMAKE) -S . -B $(Q15_BUILD) -DCMAKE_BUILD_TYPE=Release -DXOMBE_WERROR=$(WERROR) -DXOMBE_KERNELS=q15 $(CMAKE_FLAGS)
	$(CMAKE) --build $(Q15_BUILD) -j$(JOBS)
	cd $(Q15_BUILD) && ctest --output-on-failure -j$(JOBS)

conformance: all
	sh test/conformance.sh $(BUILD) $(PYTHON)

conformance-dstar: all
	sh test/conformance_dstar.sh $(BUILD) $(PYTHON)

clean:
	rm -rf $(BUILD) $(SAN_BUILD) $(Q15_BUILD)
