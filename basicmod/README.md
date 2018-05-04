# 1. Description 

The _basicmod_ library contains Fortran modules to ease several programming tasks like: error and warnings reports, work with sets, type conversions, array allocation, shell calls, command argument management, XML and VTK input/output, ...

# 2. Installation

The prerequisites are: 
- The program _make_ (for Windows it can be found as _mingw32-make_ in the MinGW distribution).
- A Fortran 2003 compiler. At the present time, only the GNU Fortran compiler, _gfortran_, and the Intel Fortran compiler, _ifort_, are supported.

After install the prerequisites:
- Go to the [basicmod webpage](https://sourceforge.net/projects/basicmod/).
- Download the compressed file and uncompress it in the installation folder.
- Open a terminal in Linux or Mac OS X, or a Command Window in Windows, go to the installation folder and type:
```shell
  make -f Makefile.<distribution>.<compiler>.<os>
```
where 
 - _&lt;distribution&gt;_ can be _static_ or _dynamic_.
 - _&lt;compiler&gt;_  can be _gfortran_ or _ifort_.
 - _&lt;os&gt;_ can be _linux_, _windows_ or _osx_. For Mac OS X, some Makefiles are provided, indicating which version they were tested for. Be aware that in Windows, MinGW distribution can use _mingw32-make_ instead of _make_.

The library and the .mod files are automatically moved to folders _lib/_  and _include/_, respectively.

When compiling for several compilers, clean the object files:
```shell
  make -f Makefile.<distribution>.<compiler>.<os> clean
```

To delete the library, the .mod and the object files, execute:
```shell
  make -f Makefile.<distribution>.<compiler>.<os> cleanlib
```

# 3. Usage

Detailed help about thye _basicmod_ library can be consulted in _./doc/index.html_.

In order to test the library, follow these steps:
- Open a terminal in Linux or Mac OS X, or a Command Window in Windows, go to the folder _test/_ and type:
```shell
  make -f Makefile.<distribution>.<compiler>.<os>
```
where _&lt;distribution&gt;_, _&lt;compiler&gt;_ and _&lt;os&gt;_ have the meaning given in the section above.

- When linking with dynamic libraries in Linux or Mac OS X, add the _lib/_ path to the variable LD\_LIBRARY\_PATH:
```shell 
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path-to-lib>
```
- When linking with dynamic libraries in Windows and _gfortran_, add the _lib/_ path to to environmental variable LD_LIBRARY_PATH. Specific instructions apply for _ifort_.
- Finally, execute the _test_ program.

# 4. Testing

This code has been tested with the following configurations:
 - Ubuntu  14.04, gfortran 4.8.4,         make 3.81 
 - CentOS   6.5,  ifort   14.0.3,         make 3.81
 - Windows 10,    gfortran 5.3.0 (MinGW), make 3.82.90
 - Windows 10,    ifort   17.0.2,         make 3.82.90

Errors have been reported with _gfortran_ version 4.6.1 and below.

# 5. License

Copyright (C) 2010-2018 Universidade de Santiago de Compostela

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.


