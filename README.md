# 1. Description

Program feconv converts finite element (FE) mesh files between several formats; it can also transform the FE type of the mesh and/or perform a bandwidth optimization. Some conversion capabilities are also present for mesh fields. Please, visit the EXAMPLES section, in the help invoked by feconv -h, to know more details.
Back to top

# 2. Installation

To install this program, you must have previously installed make in your system and a Fortran 2003 compiler. At the present time, only the GNU Fortran compiler, gfortran, and the Intel Fortran compiler, ifort, are supported.
Go to http://sourceforge.net/projects/feconv/ and download the package feconv_<date>.tar.gz, where <date> is the date of the release.
Open a terminal in Linux or Mac OS X, or a Command Window in Windows, go to the installation folder and type:

        make -f Makefile.<compiler>.<os>

where <compiler> can be "gfortran" or "ifort" and <os> can be "linux" or "windows" (for Mac OS X, "linux" is the valid option).


# 3. Supported formats
## The available input mesh formats are:

    ANSYS (.msh)
    I-Deas Universal (.unv)
    VTK-XML Unstructured Grid (.vtu)
    MD Nastran input file (.bdf)
    COMSOL mesh file (.mphtxt)
    FLUX mesh file (.pf3)
    Modulef-like Formatted Mesh (.mfm)
    Modulef-like Unformatted Mesh (.mum)
    FreFem++ Tetrahedral and/or Triangular Lagrange P1 Mesh (.msh)
    FreFem++ Tetrahedral Lagrange P1 Mesh (.mesh)
    Gmsh MSH ASCII (.msh)

## The available output mesh formats are:

    ANSYS (.msh)
    I-Deas Universal (.unv)
    VTK-XML Unstructured Grid (.vtu)
    COMSOL mesh file (.mphtxt)
    FLUX mesh file (.pf3)
    Modulef-like Formatted Mesh (.mfm)
    Modulef-like Unformatted Mesh (.mum)
    FreFem++ Tetrahedral and/or Triangular Lagrange P1 Mesh (.msh)
    FreFem++ Tetrahedral Lagrange P1 Mesh (.mesh)
    Gmsh MSH ASCII (.msh)

## The available field formats are:

    I-Deas Universal (.unv)
    FLUX mesh file (.pf3)
    FLUX field file (.dex)
    Modulef-like Formatted Field (.mff)
    Modulef-like Unformatted Field (.muf)
    ANSYS interpolation file (.ip)

# 4. License

«Copyright 2012 Iban Constenla, Victor Sande, Francisco Pena»

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

