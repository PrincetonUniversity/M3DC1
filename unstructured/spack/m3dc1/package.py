# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class M3dc1(CMakePackage):
    """M3D-C1"""

    homepage = "https://www.scorec.rpi.edu/pumi"
    git = "https://github.com/PrincetonUniversity/M3DC1.git"

    maintainers("changliu777")

    # We will use the scorec/core master branch as the 'nightly' version
    # of pumi in spack.  The master branch is more stable than the
    # scorec/core develop branch and we prefer not to expose spack users
    # to the added instability.
    version("master", submodules=True, branch="master")
    # version(
    #     "2.2.8", submodules=True, commit="736bb87ccd8db51fc499a1b91e53717a88841b1f"
    # )  # tag 2.2.8

    variant("3d", default=False, description="Enable 3D")
    variant("complex", default=False, description="Enable complex")
    # variant("zoltan", default=False, description="Enable Zoltan Features")
    # variant("fortran", default=False, description="Enable FORTRAN interface")
    # variant("testing", default=False, description="Enable all tests")
    # variant(
    #     "simmodsuite",
    #     default="none",
    #     values=("none", "base", "kernels", "full"),
    #     description="Enable Simmetrix SimModSuite Support: 'base' enables "
    #     "the minimum set of functionality, 'kernels' adds CAD kernel "
    #     "support to 'base', and 'full' enables all functionality.",
    # )
    # variant(
    #     "simmodsuite_version_check",
    #     default=True,
    #     description="Enable check of Simmetrix SimModSuite version. "
    #     "Disable the check for testing new versions.",
    # )

    conflicts("+complex", when="+3d")

    depends_on("pumi")
    depends_on("cmake@3:", type="build")
    depends_on("zoltan")
    depends_on("petsc~complex", when="~complex")
    depends_on("petsc+complex", when="+complex")
    depends_on("hdf5")
    depends_on("gsl")
    depends_on("fftw")

    def cmake_args(self):
        spec = self.spec

        args = [
            "-DCMAKE_C_COMPILER=%s" % spec["mpi"].mpicc,
            "-DCMAKE_CXX_COMPILER=%s" % spec["mpi"].mpicxx,
            "-DCMAKE_Fortran_COMPILER=%s" % spec["mpi"].mpifc,
            "-DENABLE_PETSC=ON",
            "-DENABLE_ZOLTAN=ON",
            self.define_from_variant("ENABLE_3D", "3d"),
            self.define_from_variant("ENABLE_COMPLEX", "complex"),
        ]
        return args
