# fastcpd C++ example

This example uses the standalone C++ header API and the shared
`tests/fixtures/manifest.tsv` contract:

```cpp
#include <fastcpd/fastcpd.h>
```

Build from this directory with plain Make:

```sh
make
./mean_change
```

The example selects the mean detector row from the shared manifest and loads
its CSV relative to the manifest. When run from the repository root (as CI
does), it discovers `tests/fixtures/manifest.tsv` automatically. When run from
this directory, it searches `../../tests/fixtures/manifest.tsv`; an explicit
fixture directory or manifest path may be supplied as the first argument:

```sh
./mean_change ../../tests/fixtures
```

If Abseil is not installed in a standard include path, pass its include
directory explicitly:

```sh
make ABSL_INCLUDE=/path/to/abseil/include
```

If Armadillo is not installed in a standard include path, pass its include
directory explicitly:

```sh
make ARMADILLO_INCLUDE=/path/to/armadillo/include
```

The Makefile also auto-detects the R `abseil` and `RcppArmadillo` package
include directories when R is available. For CMake consumers, configure the
repository root and link the interface target:

```cmake
find_package(fastcpd REQUIRED)
target_link_libraries(my_target PRIVATE fastcpd::fastcpd)
```

The CMake package requires Abseil `20260526` or newer. If Abseil is installed
under a non-standard prefix, add that prefix to `CMAKE_PREFIX_PATH` or set
`absl_DIR` to the directory containing `abslConfig.cmake`.
