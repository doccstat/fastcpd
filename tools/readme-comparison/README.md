# README comparison reproduction

Run from the canonical fastcpd source directory. R needs fastcpd, mosum,
changepoint, fpop, callr, ggplot2, knitr, and rmarkdown.

## Dependencies and native build

Tested on Linux ARM64 with GCC 13.3, fastcpd 1.3.0, Python 3.12,
NumPy 2.5.2, and Numba 0.67.0.

```sh
python -m pip install 'ruptures==1.1.10' 'skchange[numba]==0.18.0' 'sdt-python==20.1.4'
```

Build and install fastcpd in Release mode with Armadillo and Abseil.
The benchmark compiles fpop's unmodified C++ core without its R wrapper;
its LGPL-2.1-or-later sources stay in a separate checkout.

```sh
git clone https://github.com/cran/fpop.git /tmp/fastcpd-readme-fpop
git -C /tmp/fastcpd-readme-fpop switch --detach 855c14826b568665f8318e270b49ea6c267f8a31
cmake -S tools/readme-comparison -B /tmp/fastcpd-readme-build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH='/path/to/fastcpd-install;/path/to/abseil-install' \
  -DFPOP_SOURCE_DIR=/tmp/fastcpd-readme-fpop
cmake --build /tmp/fastcpd-readme-build --parallel 2
```

## Render

Set `R_LIBS` if fastcpd is installed in a separate R library.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
FASTCPD_README_PYTHON=/path/to/venv/bin/python \
FASTCPD_README_CPP=/tmp/fastcpd-readme-build/readme_comparison \
Rscript -e 'rmarkdown::render("README.Rmd")'
```

The render regenerates all three figures and passes a million-observation
input to the native executable through a temporary file. Benchmark settings and
timing boundaries are defined in `README.Rmd` and `comparison.cc`.
