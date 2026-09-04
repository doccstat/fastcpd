"""Tests for source, artifact, and tag version consistency."""

from email.message import Message
import importlib.util
from io import BytesIO
from pathlib import Path
import sys
import tarfile
import zipfile

import pytest


CHECKER_PATH = (
    Path(__file__).resolve().parents[1] / "python" / "_check_versions.py"
)
if not CHECKER_PATH.exists():
    pytest.skip(
        "version checker is a source-distribution/repository tool",
        allow_module_level=True,
    )
SPEC = importlib.util.spec_from_file_location("fastcpd_version_check", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
version_check = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = version_check
SPEC.loader.exec_module(version_check)


def _metadata(version):
    message = Message()
    message["Metadata-Version"] = "2.4"
    message["Name"] = "fastcpd"
    message["Version"] = version
    return message.as_bytes()


def _write_sources(root, version="1.3.0"):
    (root / "python").mkdir()
    (root / "python" / "__init__.py").write_text(
        f'__version__ = "{version}"\n', encoding="utf-8"
    )
    (root / "DESCRIPTION").write_text(
        f"Package: fastcpd\nVersion: {version}\n", encoding="utf-8"
    )
    (root / "CMakeLists.txt").write_text(
        f"project(fastcpd VERSION {version} LANGUAGES CXX)\n",
        encoding="utf-8",
    )


def _write_wheel(artifact_dir, version="1.3.0"):
    path = artifact_dir / f"fastcpd-{version}-py3-none-any.whl"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"fastcpd-{version}.dist-info/METADATA", _metadata(version))
    return path


def _write_sdist(artifact_dir, version="1.3.0"):
    path = artifact_dir / f"fastcpd-{version}.tar.gz"
    contents = _metadata(version)
    info = tarfile.TarInfo(f"fastcpd-{version}/PKG-INFO")
    info.size = len(contents)
    with tarfile.open(path, "w:gz") as archive:
        archive.addfile(info, BytesIO(contents))
    return path


def test_unified_sources_artifacts_and_tag_share_one_version(tmp_path):
    _write_sources(tmp_path)
    artifact_dir = tmp_path / "dist"
    artifact_dir.mkdir()
    _write_wheel(artifact_dir)
    _write_sdist(artifact_dir)

    observations = version_check.check_versions(
        tmp_path,
        artifact_dirs=(artifact_dir,),
        tag="v1.3.0",
        unified=True,
        require_wheel=True,
        require_sdist=True,
    )
    assert observations
    assert {item.version for item in observations} == {"1.3.0"}


def test_unified_source_mismatch_fails_with_the_source_name(tmp_path):
    _write_sources(tmp_path)
    (tmp_path / "CMakeLists.txt").write_text(
        "project(fastcpd VERSION 1.2.9 LANGUAGES CXX)\n",
        encoding="utf-8",
    )

    with pytest.raises(version_check.VersionCheckError, match="CMakeLists.txt"):
        version_check.check_versions(tmp_path, unified=True)


def test_artifact_metadata_and_release_tag_mismatches_fail(tmp_path):
    _write_sources(tmp_path)
    artifact_dir = tmp_path / "dist"
    artifact_dir.mkdir()
    wheel = _write_wheel(artifact_dir)
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("fastcpd-1.3.0.dist-info/METADATA", _metadata("1.3.1"))

    with pytest.raises(version_check.VersionCheckError, match="wheel metadata"):
        version_check.check_versions(tmp_path, artifact_dirs=(artifact_dir,))
    with pytest.raises(version_check.VersionCheckError, match="release tag"):
        version_check.check_versions(tmp_path, tag="v1.3.1")
