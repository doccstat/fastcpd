load("//bazel:repositories.bzl", "com_github_pybind11_pybind")

com_github_pybind11_pybind()

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "rules_python",
    sha256 = "9f9f3b300a9264e4c77999312ce663be5dee9a56e361a1f6fe7ec60e1beef9a3",
    strip_prefix = "rules_python-1.4.1",
    url = "https://github.com/bazel-contrib/rules_python/releases/download/1.4.1/rules_python-1.4.1.tar.gz",
)

load("@rules_python//python:repositories.bzl", "py_repositories", "python_register_toolchains")

py_repositories()

python_register_toolchains(
    name = "python_3_10",
    python_version = "3.10",
)

load("@rules_python//python:pip.bzl", "pip_parse")

pip_parse(
    # (Optional) You can set an environment in the pip process to control its
    # behavior. Note that pip is run in "isolated" mode so no PIP_<VAR>_<NAME>
    # style env vars are read, but env vars that control requests and urllib3
    # can be passed
    # environment = {"HTTPS_PROXY": "http://my.proxy.fun/"},
    name = "pypi",

    # Requirement groups allow Bazel to tolerate PyPi cycles by putting dependencies
    # which are known to form cycles into groups together.
    # experimental_requirement_cycles = {
    #     "sphinx": [
    #         "sphinx",
    #         "sphinxcontrib-qthelp",
    #         "sphinxcontrib-htmlhelp",
    #         "sphinxcontrib-devhelp",
    #         "sphinxcontrib-applehelp",
    #         "sphinxcontrib-serializinghtml",
    #     ],
    # },
    # (Optional) You can provide extra parameters to pip.
    # Here, make pip output verbose (this is usable with `quiet = False`).
    # extra_pip_args = ["-v"],

    # (Optional) You can exclude custom elements in the data section of the generated BUILD files for pip packages.
    # Exclude directories with spaces in their names in this example (avoids build errors if there are such directories).
    #pip_data_exclude = ["**/* */**"],

    # (Optional) You can provide a python_interpreter (path) or a python_interpreter_target (a Bazel target, that
    # acts as an executable). The latter can be anything that could be used as Python interpreter. E.g.:
    # 1. Python interpreter that you compile in the build file (as above in @python_interpreter).
    # 2. Pre-compiled python interpreter included with http_archive
    # 3. Wrapper script, like in the autodetecting python toolchain.
    #
    # Here, we use the interpreter constant that resolves to the host interpreter from the default Python toolchain.
    python_interpreter_target = "@python_3_10_host//:python",

    # (Optional) You can set quiet to False if you want to see pip output.
    #quiet = False,
    requirements_lock = "//:requirements_lock.txt",
)

load("@pypi//:requirements.bzl", "install_deps")

# Initialize repositories for all packages in requirements_lock.txt.
install_deps()

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "rules_java",
    sha256 = "1558508fc6c348d7f99477bd21681e5746936f15f0436b5f4233e30832a590f9",
    urls = [
        "https://github.com/bazelbuild/rules_java/releases/download/8.12.0/rules_java-8.12.0.tar.gz",
    ],
)

load("@rules_java//java:rules_java_deps.bzl", "rules_java_dependencies")

rules_java_dependencies()

load("@bazel_features//:deps.bzl", "bazel_features_deps")

bazel_features_deps()

# note that the following line is what is minimally required from protobuf for the java rules
# consider using the protobuf_deps() public API from @com_google_protobuf//:protobuf_deps.bzl
load("@com_google_protobuf//bazel/private:proto_bazel_features.bzl", "proto_bazel_features")  # buildifier: disable=bzl-visibility

proto_bazel_features(name = "proto_bazel_features")

# register toolchains
load("@rules_java//java:repositories.bzl", "rules_java_toolchains")

rules_java_toolchains()

load("//bazel/python:python_native.bzl", "local_python_deps")

local_python_deps()

http_archive(
    name = "armadillo_headers",
    build_file_content = """
cc_library(
    name = "armadillo_header",
    hdrs = glob(["include/armadillo", "include/armadillo_bits/*.hpp"]),
    includes = ["include/"],
    visibility = ["//visibility:public"],
)
""",
    sha256 = "023242fd59071d98c75fb015fd3293c921132dc39bf46d221d4b059aae8d79f4",
    strip_prefix = "armadillo-14.4.0",
    urls = ["https://sourceforge.net/projects/arma/files/armadillo-14.4.0.tar.xz"],
)

http_archive(
    name = "abseil-cpp",
    sha256 = "7262daa7c1711406248c10f41026d685e88223bc92817d16fb93c19adb57f669",
    strip_prefix = "abseil-cpp-20250512.0",
    urls = ["https://github.com/abseil/abseil-cpp/releases/download/20250512.0/abseil-cpp-20250512.0.tar.gz"],
)

# Pre-built OpenBLAS for Windows x64. Provides libopenblas.lib (import lib)
# and libopenblas.dll for BLAS/LAPACK on MSVC builds.
http_archive(
    name = "openblas_windows",
    build_file_content = """
cc_library(
    name = "openblas",
    hdrs = glob(["include/*.h"]),
    includes = ["include"],
    srcs = ["lib/libopenblas.lib"],
    data = ["bin/libopenblas.dll"],
    visibility = ["//visibility:public"],
)
""",
    sha256 = "b42a74d1c9c77bdab2cf2688031b9bc4a322ade71c549427f4950d85fd590fca",
    url = "https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29_x64.zip",
)

http_archive(
    name = "com_google_protobuf",
    sha256 = "1451b03faec83aed17cdc71671d1bbdfd72e54086b827f5f6fd02bf7a4041b68",
    strip_prefix = "protobuf-30.1",
    url = "https://github.com/protocolbuffers/protobuf/releases/download/v30.1/protobuf-30.1.tar.gz",
)

load("@com_google_protobuf//:protobuf_deps.bzl", "protobuf_deps")

protobuf_deps()

http_archive(
    name = "googleapis",
    sha256 = "e1fd1bac378a21a2557439dff96114f8bb24e985630de455116a33c7c60aa619",
    strip_prefix = "googleapis-111b7383752255d1849a8d3b7259ed735acc4f97",
    url = "https://github.com/googleapis/googleapis/archive/111b7383752255d1849a8d3b7259ed735acc4f97.tar.gz",
)

http_archive(
    name = "googletest",
    sha256 = "78c676fc63881529bf97bf9d45948d905a66833fbfa5318ea2cd7478cb98f399",
    strip_prefix = "googletest-1.16.0",
    url = "https://github.com/google/googletest/releases/download/v1.16.0/googletest-1.16.0.tar.gz",
)
