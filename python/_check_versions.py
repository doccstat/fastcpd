#!/usr/bin/env python3
"""Verify fastcpd source, artifact, and release-tag versions.

This private repository tool intentionally uses only the Python standard
library so the source version can be checked before the fastcpd extension or
build dependencies are installed.
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
from email.parser import BytesParser
from pathlib import Path
import re
import sys
import tarfile
import zipfile


_DESCRIPTION_VERSION = re.compile(r"^Version:\s*(\S+)\s*$", re.MULTILINE)
_CMAKE_PROJECT = re.compile(
    r"project\s*\(\s*fastcpd\b(?P<body>.*?)\)",
    re.IGNORECASE | re.DOTALL,
)
_CMAKE_VERSION = re.compile(r"\bVERSION\s+([^\s)]+)", re.IGNORECASE)


class VersionCheckError(RuntimeError):
    """Raised when version-bearing files do not describe one release."""


@dataclass(frozen=True)
class VersionObservation:
    """One version read from a source, artifact, or tag."""

    label: str
    version: str


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError as error:
        raise VersionCheckError(f"cannot read {path}: {error}") from error


def read_python_version(path: Path) -> str:
    """Read ``__version__`` without importing the compiled package."""
    try:
        module = ast.parse(_read_text(path), filename=str(path))
    except SyntaxError as error:
        raise VersionCheckError(f"cannot parse {path}: {error}") from error
    for node in module.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        if not any(
            isinstance(target, ast.Name) and target.id == "__version__"
            for target in targets
        ):
            continue
        value = node.value
        if isinstance(value, ast.Constant) and isinstance(value.value, str):
            return value.value
        raise VersionCheckError(f"{path}: __version__ must be a string literal")
    raise VersionCheckError(f"{path}: no __version__ assignment found")


def read_description_version(path: Path) -> str:
    """Read the R package version from DESCRIPTION."""
    match = _DESCRIPTION_VERSION.search(_read_text(path))
    if match is None:
        raise VersionCheckError(f"{path}: no Version field found")
    return match.group(1)


def read_cmake_version(path: Path) -> str:
    """Read the standalone C++ project version from CMakeLists.txt."""
    project = _CMAKE_PROJECT.search(_read_text(path))
    if project is None:
        raise VersionCheckError(f"{path}: no fastcpd project() call found")
    version = _CMAKE_VERSION.search(project.group("body"))
    if version is None:
        raise VersionCheckError(f"{path}: fastcpd project() has no VERSION")
    return version.group(1)


def _metadata_version(contents: bytes, source: str) -> str:
    metadata = BytesParser().parsebytes(contents)
    if metadata.get("Name", "").lower().replace("_", "-") != "fastcpd":
        raise VersionCheckError(f"{source}: package Name is not fastcpd")
    version = metadata.get("Version")
    if not version:
        raise VersionCheckError(f"{source}: package Version is missing")
    return version


def read_wheel_versions(path: Path) -> list[VersionObservation]:
    """Read the version from a wheel filename and its METADATA."""
    parts = path.name.removesuffix(".whl").split("-")
    if len(parts) < 5 or parts[0].lower().replace("_", "-") != "fastcpd":
        raise VersionCheckError(f"{path}: malformed fastcpd wheel filename")
    filename_version = parts[1]
    try:
        with zipfile.ZipFile(path) as archive:
            metadata_names = [
                name for name in archive.namelist()
                if name.endswith(".dist-info/METADATA")
            ]
            if len(metadata_names) != 1:
                raise VersionCheckError(
                    f"{path}: expected one .dist-info/METADATA entry, "
                    f"found {len(metadata_names)}"
                )
            metadata_version = _metadata_version(
                archive.read(metadata_names[0]), f"{path}:{metadata_names[0]}"
            )
    except (OSError, zipfile.BadZipFile) as error:
        raise VersionCheckError(f"cannot inspect wheel {path}: {error}") from error
    return [
        VersionObservation(f"wheel filename {path.name}", filename_version),
        VersionObservation(f"wheel metadata {path.name}", metadata_version),
    ]


def read_sdist_versions(path: Path) -> list[VersionObservation]:
    """Read the version from an sdist filename and its PKG-INFO."""
    prefix = "fastcpd-"
    suffix = ".tar.gz"
    if not path.name.startswith(prefix) or not path.name.endswith(suffix):
        raise VersionCheckError(f"{path}: malformed fastcpd sdist filename")
    filename_version = path.name[len(prefix):-len(suffix)]
    try:
        with tarfile.open(path, mode="r:gz") as archive:
            metadata_members = [
                member for member in archive.getmembers()
                if Path(member.name).name == "PKG-INFO" and member.isfile()
            ]
            if len(metadata_members) != 1:
                raise VersionCheckError(
                    f"{path}: expected one PKG-INFO entry, "
                    f"found {len(metadata_members)}"
                )
            extracted = archive.extractfile(metadata_members[0])
            if extracted is None:
                raise VersionCheckError(f"{path}: cannot read PKG-INFO")
            metadata_version = _metadata_version(
                extracted.read(), f"{path}:{metadata_members[0].name}"
            )
    except (OSError, tarfile.TarError) as error:
        raise VersionCheckError(f"cannot inspect sdist {path}: {error}") from error
    return [
        VersionObservation(f"sdist filename {path.name}", filename_version),
        VersionObservation(f"sdist metadata {path.name}", metadata_version),
    ]


def check_versions(
    root: Path,
    artifact_dirs: tuple[Path, ...] = (),
    tag: str | None = None,
    unified: bool = False,
    require_wheel: bool = False,
    require_sdist: bool = False,
) -> list[VersionObservation]:
    """Return all observations after verifying they match Python's version."""
    root = root.resolve()
    python_path = root / "python" / "__init__.py"
    expected = read_python_version(python_path)
    observations = [VersionObservation(str(python_path), expected)]

    if unified:
        observations.extend([
            VersionObservation(
                str(root / "DESCRIPTION"),
                read_description_version(root / "DESCRIPTION"),
            ),
            VersionObservation(
                str(root / "CMakeLists.txt"),
                read_cmake_version(root / "CMakeLists.txt"),
            ),
        ])

    wheels: list[Path] = []
    sdists: list[Path] = []
    for artifact_dir in artifact_dirs:
        artifact_dir = artifact_dir.resolve()
        wheels.extend(sorted(artifact_dir.rglob("fastcpd-*.whl")))
        sdists.extend(sorted(artifact_dir.rglob("fastcpd-*.tar.gz")))
    if require_wheel and not wheels:
        raise VersionCheckError("no fastcpd wheel found in the artifact paths")
    if require_sdist and not sdists:
        raise VersionCheckError("no fastcpd sdist found in the artifact paths")
    for wheel in wheels:
        observations.extend(read_wheel_versions(wheel))
    for sdist in sdists:
        observations.extend(read_sdist_versions(sdist))

    if tag is not None:
        expected_tag = f"v{expected}"
        if tag != expected_tag:
            raise VersionCheckError(
                f"release tag {tag!r} does not match {expected_tag!r}"
            )
        observations.append(VersionObservation("release tag", tag[1:]))

    mismatches = [item for item in observations if item.version != expected]
    if mismatches:
        details = "; ".join(
            f"{item.label} reports {item.version!r}" for item in mismatches
        )
        raise VersionCheckError(
            f"expected every version to be {expected!r}; {details}"
        )
    return observations


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="canonical source root (default: parent of python/)",
    )
    parser.add_argument(
        "--artifacts",
        type=Path,
        action="append",
        default=[],
        help="directory searched recursively for wheels and sdists",
    )
    parser.add_argument("--tag", help="release tag, including the leading v")
    parser.add_argument(
        "--unified",
        action="store_true",
        help="also require R DESCRIPTION and CMake project versions to match",
    )
    parser.add_argument("--require-wheel", action="store_true")
    parser.add_argument("--require-sdist", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        observations = check_versions(
            root=args.root,
            artifact_dirs=tuple(args.artifacts),
            tag=args.tag,
            unified=args.unified,
            require_wheel=args.require_wheel,
            require_sdist=args.require_sdist,
        )
    except VersionCheckError as error:
        print(f"version consistency check failed: {error}", file=sys.stderr)
        return 1
    for item in observations:
        print(f"{item.label}: {item.version}")
    print(f"version consistency check passed for {observations[0].version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
