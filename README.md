# 1. Description

FEconv can convert finite element (FE) mesh written in several commercial file formats. It can also transform the FE type and/or perform some bandwidth optimizations. Some conversion capabilities are also present for mesh fields. Please type `feconv -h` and see EXAMPLES section to know more details, or visit <a href="http://victorsndvg.github.io/FEconv/">http://victorsndvg.github.io/FEconv</a> for more information.

# 2. Installation

## As a standalone program

The prerequisites are: 
 - the program _make_ (for Windows it can be found as _mingw32-make_ in the MinGW distribution) and 
 - a Fortran 2003 compiler; at the present time, only the GNU Fortran compiler, _gfortran_, and the Intel Fortran compiler, _ifort_, are supported in FEconv.
 
After install the prerequisites, go to https://github.com/victorsndvg/FEconv, download the ZIP file and unzip it in the installation folder; open a terminal in Linux or Mac OS X, or a Command Window in Windows, go to the installation folder and type:
```shell
  make -f Makefile.<compiler>.<os>
```
where _\<compiler\>_  can be _gfortran_ or _ifort_ and _\<os\>_ can be _linux_, _windows_ or _osx_. For Mac OS X, some Makefiles are provided, indicating for which version they were tested. Be aware that in Windows, MinGW distribution uses _mingw32-make_ instead of _make_.

## As a library

Prebuilt libraries and header files are located in folders _lib/_  and _include/_. 

If you want to build the libraries by your own, you must install the prerequisites previously mentioned. Then open a terminal in Linux or Mac Mac OS X, or a Command Window in Windows, go to the installation folder and type:
```shell
  make -f Makefile.makelib.<compiler>.<os>
```
where _\<compiler\>_ can be _gfortran_ or _ifort_ and _\<os\>_ can be _linux_, _windows_ or _osx_. For Mac OS X, some Makefiles are provided, indicating for which version they were tested. Be aware that in Windows, MinGW distribution uses _mingw32-make_ instead of _make_.

# 3. Usage

## As a standalone program

Please execute `feconv -h` to see the command line options and some examples of use, or visit the [FEconv help webpage](http://victorsndvg.github.io/FEconv/). 

## As a library

Please inspect the folder _testlib/_ to see an example of library use. 

# 4. Supported formats
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
    VTK-XML Unstructured Grid (.vtu)
    FLUX mesh file (.pf3)
    FLUX field file (.dex)
    Modulef-like Formatted Field (.mff)
    Modulef-like Unformatted Field (.muf)
    ANSYS interpolation file (.ip)

# 5. License

Copyright (C) 2010-2018 Universidade de Santiago de Compostela

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

