# 1. Description

FEconv can convert finite element (FE) mesh written in several commercial file formats. It can also transform the FE type and/or perform some bandwidth optimizations. Some conversion capabilities are also present for mesh fields. Please type `feconv -h` and see EXAMPLES section to know more details, or visit the [FEconv help webpage](http://victorsndvg.github.io/FEconv) for more information.

# 2. Installation

## As a standalone program

The prerequisites are: 
- The program _make_ (for Windows it can be found as _mingw32-make_ in the MinGW distribution).
- A Fortran 2003 compiler. At the present time, only the GNU Fortran compiler, _gfortran_, and the Intel Fortran compiler, _ifort_, are supported.
 
After install the prerequisites:
- Go to the [FEconv webpage](https://github.com/victorsndvg/FEconv).
- Download the compressed file and uncompress it in the installation folder.
- Open a terminal in Linux or Mac OS X, or a Command Window in Windows, go to the installation folder and type:
```shell
  make -f Makefile.<compiler>.<os>
```
where 
 - _&lt;compiler&gt;_  can be _gfortran_ or _ifort_.
 - _&lt;os&gt;_ can be _linux_, _windows_ or _osx_. For Mac OS X, some Makefiles are provided, indicating which version they were tested for. Be aware that in Windows, MinGW distribution can use _mingw32-make_ instead of _make_.

To delete the FEconv executable, its .mod and object files, execute:
```shell
  make -f Makefile.<compiler>.<os> clean
```

To delete the basicmod library, its .mod and object files, execute:
```shell
  make -f Makefile.<compiler>.<os> cleanlib
```

## As a library

- The command to execute in the terminal or Command Window is slightly different from the previous one:
```shell
  make -f Makefile.makelib.<distribution>.<compiler>.<os>
```
where _&lt;distribution&gt;_ can be _static_ or _dynamic_.

The library and the .mod files are automatically moved to folders _lib/_  and _include/_, respectively.

When compiling for several compilers, clean the object files:
```shell
  make -f Makefile.makelib.<distribution>.<compiler>.<os> clean
```

To delete the FEconv (and basicmod) libraries, the .mod and the object files, execute:
```shell
  make -f Makefile.makelib.<distribution>.<compiler>.<os> cleanlib
```
# 3. Usage

## As a standalone program

Please execute `feconv -h` to see the command line options and some examples of use, or visit the [FEconv help webpage](http://victorsndvg.github.io/FEconv). 

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

Copyright (C) 2010-2020 Universidade de Santiago de Compostela

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

