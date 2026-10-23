from __future__ import annotations

from pathlib import Path
import shutil

from setuptools import Extension, find_namespace_packages, setup
from setuptools.command.build_ext import build_ext

PREBUILT = Path(__file__).parent / "prebuilt"

NANOBIND_MODULES = [
    "lillymol",
    "lillymol_bdb",
    "lillymol_gfp_server",
]

COMPATIBILITY_MODULES = [
    "lillymol_atom",
    "lillymol_bond",
    "lillymol_fingerprint",
    "lillymol_io",
    "lillymol_query",
    "lillymol_reaction",
    "lillymol_ring",
    "lillymol_set_of_atoms",
    "lillymol_standardise",
    "lillymol_tools",
    "lillymol_tsubstructure",
]


class PrebuiltExtension(Extension):
  """A setuptools Extension whose binary is built by Bazel."""

  def __init__(self, name: str):
    super().__init__(name, sources=[])


class BuildPrebuiltExtensions(build_ext):
  """Copy staged Bazel-built extension modules into the wheel build tree."""

  def run(self):
    if not PREBUILT.is_dir():
      raise RuntimeError(
          f"{PREBUILT} does not exist. Run scripts/stage_wheel_files.sh first.")

    super().run()
    self._copy_private_shared_libraries()

  def build_extension(self, ext):
    source = PREBUILT / f"{ext.name}.so"
    if not source.is_file():
      raise RuntimeError(f"Missing prebuilt extension module {source}")

    destination = Path(self.get_ext_fullpath(ext.name))
    destination.parent.mkdir(parents=True, exist_ok=True)
    self.copy_file(str(source), str(destination))

  def _copy_private_shared_libraries(self):
    build_lib = Path(self.build_lib)
    build_lib.mkdir(parents=True, exist_ok=True)
    for source in sorted(PREBUILT.glob("lib*.so*")):
      destination = build_lib / source.name
      self.copy_file(str(source), str(destination))


setup(
    packages=find_namespace_packages(where="prebuilt"),
    package_dir={"": "prebuilt"},
    py_modules=COMPATIBILITY_MODULES,
    ext_modules=[PrebuiltExtension(name) for name in NANOBIND_MODULES],
    cmdclass={"build_ext": BuildPrebuiltExtensions},
    zip_safe=False,
)
