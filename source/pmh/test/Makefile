#################################################################
# Makefile created using the tool 'Creamake'
# 
# Creamake is distributed under the GNU GPL license
# Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
# Download page: http://sourceforge.net/projects/creamake/
#################################################################
 
#################################################################
# User-modifiable options
#################################################################
 
# SOURCE AND COMMONS FOLDERS (separated by spaces)
dir_fuentes = . .. ../../basicmod ../../basicmod/gfortran ../../basicmod/args \
../../basicmod/alloc ../../basicmod/vtk ../../mfm ../../vtu
 
# OBJECT AND .MOD FOLDER
dir_objetos = object
 
# MAIN SOURCE FILE (include relative path from folder where Makefile is)
condir_principal = ./test_pmh.f90
 
# EXECUTABLE NAME 
ejecutable = test_pmh.exe
 
# NEEDED TO convert ejecutable THE DEFAULT RULE: 
$(ejecutable): $(condir_principal) 
 
# MODULES
modulos = module_fe_database_pmh.f90 module_compiler_gfortran.f90 \
module_os_dependant.f90 module_report.f90 module_set.f90 module_convers.f90 \
module_files.f90 module_feed.f90 module_args.f90 module_alloc_char_r1.f90 \
module_alloc_int_r1.f90 module_alloc_int_r2.f90 module_alloc_log_r1.f90 \
module_alloc_real64_r1.f90 module_alloc_real64_r2.f90 module_alloc.f90 \
module_pmh.f90 module_system.f90 IR_Precision.f90 Lib_VTK_IO.f90 \
LIB_VTK_IO_READ.f90 module_writevtu.f90 module_mfm.f90 module_vtu.f90
 
# MODULE DEPENDENCIES
# if pru1 depends on pru2... pru1.o: pru2.o
module_os_dependant.o: module_compiler_gfortran.o
module_report.o: module_compiler_gfortran.o module_os_dependant.o
module_set.o: module_os_dependant.o module_report.o
module_convers.o: module_compiler_gfortran.o module_os_dependant.o \
module_report.o
module_files.o: module_os_dependant.o module_report.o module_convers.o
module_feed.o: module_convers.o
module_args.o: module_compiler_gfortran.o module_os_dependant.o \
module_report.o module_convers.o
module_alloc_char_r1.o: module_os_dependant.o module_report.o
module_alloc_int_r1.o: module_os_dependant.o module_report.o
module_alloc_int_r2.o: module_os_dependant.o module_report.o \
module_alloc_int_r1.o
module_alloc_log_r1.o: module_os_dependant.o module_report.o \
module_alloc_int_r1.o
module_alloc_real64_r1.o: module_compiler_gfortran.o module_os_dependant.o \
module_report.o module_alloc_int_r1.o
module_alloc_real64_r2.o: module_compiler_gfortran.o module_os_dependant.o \
module_report.o module_alloc_real64_r1.o
module_alloc.o: module_alloc_int_r1.o module_alloc_int_r2.o \
module_alloc_real64_r1.o module_alloc_real64_r2.o module_alloc_char_r1.o \
module_alloc_log_r1.o
module_pmh.o: module_compiler_gfortran.o module_os_dependant.o module_report.o \
module_convers.o module_alloc.o module_args.o module_fe_database_pmh.o
module_system.o: module_compiler_gfortran.o module_os_dependant.o \
module_report.o module_convers.o module_alloc_char_r1.o module_files.o
Lib_VTK_IO.o: IR_Precision.o
LIB_VTK_IO_READ.o: Lib_VTK_IO.o
module_writevtu.o: Lib_VTK_IO.o
module_mfm.o: module_compiler_gfortran.o module_os_dependant.o module_report.o \
module_convers.o module_feed.o
module_vtu.o: module_compiler_gfortran.o module_report.o module_convers.o \
module_alloc_int_r1.o module_set.o module_writevtu.o
 
# INCLUDES
includes = 
 
# COMPILER
FC = gfortran
 
# COMPILER OPTIONS
FFLAGS = -J$(dir_objetos) -Wall -g -fcheck=all -fbacktrace -fall-intrinsics -L /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/lib
 
# LINKER OPTIONS
LDFLAGS = 
 
#################################################################
# Non-modifiable part
#################################################################
 
# SOURCE FOLDERS
VPATH =   $(subst ,:,$(strip $(dir_fuentes)))
vpath %.o $(dir_objetos)
 
# SOURCES
fuentes_ = $(filter %.f %.F %.for %.FOR %.f90 %.F90 %.f95 %.F95 %.f03 %.F03,$(shell ls $(dir_fuentes)))
fuentes  = $(filter-out $(notdir $(condir_principal)) $(modulos),$(fuentes_))
 
# OBJECTS
modulos_obj = $(addsuffix .o,$(basename $(modulos)))
fuentes_obj = $(addsuffix .o,$(basename $(fuentes)))
 
# OBJECTS WITH PATH
condir_modulos_obj = $(addprefix $(dir_objetos)/,$(modulos_obj))
condir_fuentes_obj = $(addprefix $(dir_objetos)/,$(fuentes_obj))
 
# COMPILATION OPTIONS
FFLAGS += $(patsubst %,-I%,$(dir_fuentes))
FFLAGS += -I$(dir_objetos)
 
# MAIN RULE
all: $(ejecutable)
 
$(ejecutable): $(includes) $(modulos_obj) $(fuentes_obj)
	$(FC) -o $(ejecutable) $(FFLAGS) $(condir_principal) $(condir_modulos_obj) $(condir_fuentes_obj) $(LDFLAGS)
 
# SOURCES RULE
$(fuentes_obj): $(includes) $(modulos_obj)
 
# RULE PATTERNS
%.o:%.f
	$(FC) -c -o $@ $(FFLAGS) $<
	@mv $@ $(dir_objetos) 
%.o:%.F
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
%.o:%.for
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
%.o:%.FOR
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
%.o:%.f90
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
%.o:%.F90
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
 
.PHONY: clean
clean:
	rm $(dir_objetos)/*.o      
	rm $(dir_objetos)/*.mod    
	rm $(ejecutable)
