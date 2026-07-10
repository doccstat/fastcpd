# fastcpd C++ example

This example uses the standalone C++ header API:

```cpp
#include <fastcpd/fastcpd.h>
```

Build from this directory with plain Make:

```sh
make
./simple
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
