#################################################################
# Makefile created using the tool 'Creamake'
# 
# Creamake is distributed under the GNU GPL license
# Author: Francisco Pena, fran.pena@usc.es
# Download page: http://sourceforge.net/projects/creamake/
#################################################################
 
#################################################################
# User-modifiable options
#################################################################
 
# SOURCE AND COMMONS FOLDERS (separated by spaces)
dir_fuentes = . ..
 
# OBJECT AND .MOD FOLDER
dir_objetos = object
 
# MAIN SOURCE FILE (include relative path from folder where Makefile is)
condir_principal = ./test_gfortran.f90
 
# EXECUTABLE NAME 
ejecutable = test
 
# NEEDED TO convert ejecutable THE DEFAULT RULE: 
$(ejecutable): $(condir_principal) 
 
# MODULES
modulos = module_compiler_gfortran.f90
 
# MODULE DEPENDENCIES
# if pru1 depends on pru2... pru1.o: pru2.o
 
# INCLUDES
includes = 
 
# COMPILER
FC = gfortran
 
# COMPILER OPTIONS
FFLAGS = -J$(dir_objetos)
 
# LINKER OPTIONS
LDFLAGS = 
 
#################################################################
# Non-modifiable part
#################################################################
 
# CARPETAS DE BUSQUEDA
VPATH =   $(subst ,:,$(strip $(dir_fuentes)))
vpath %.o $(dir_objetos)
 
# FUENTES
fuentes_ = $(filter %.f %.F %.for %.FOR %.f90 %.F90 %.f95 %.F95 %.f03 %.F03,$(shell ls $(dir_fuentes)))
fuentes  = $(filter-out $(notdir $(condir_principal)) $(modulos),$(fuentes_))
 
# OBJETOS
modulos_obj = $(addsuffix .o,$(basename $(modulos)))
fuentes_obj = $(addsuffix .o,$(basename $(fuentes)))
 
# OBJETOS EN SU CARPETA
condir_modulos_obj = $(addprefix $(dir_objetos)/,$(modulos_obj))
condir_fuentes_obj = $(addprefix $(dir_objetos)/,$(fuentes_obj))
 
# OPCIONES DE COMPILACION
FFLAGS += $(patsubst %,-I%,$(dir_fuentes))
FFLAGS += -I$(dir_objetos)
 
# REGLA PRINCIPAL
all: $(ejecutable)
 
$(ejecutable): $(includes) $(modulos_obj) $(fuentes_obj)
	$(FC) -o $(ejecutable) $(FFLAGS) $(condir_principal) $(condir_modulos_obj) $(condir_fuentes_obj) $(LDFLAGS)
 
# REGLA PARA LAS FUENTES
$(fuentes_obj): $(includes) $(modulos_obj)
 
# PATRONES DE REGLAS
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
