module basicmod
!! License: GNU GPLv2 license.
!! Author: Francisco Pena, fran.pena(at)usc.es
!! Date: 31/12/2017
!!
!! Module to call `basicmod` entities.  
!!
!! @note A program can also call other procedures with identical names than those from `basicmod`. See the following example 
!! contained in folder `test`.  
!!
!! __Example:__  
!! `module module_report_old`  
!! `implicit none`  
!! `contains`  
!! `subroutine info(msg)`  
!! `character(*), intent(in) :: msg`  
!! `print*, trim(msg)`  
!! `end subroutine`  
!! `end module`  
!!
!! `program test`  
!! `use basicmod, only: report_option, get_os, info`  
!! `use module_report_old, only: linfo => info`  
!! `implicit none`  
!! `call report_option('info', 'std')`  
!! `call info('Detected OS using basicmod: '//trim(get_os()))`  
!! `call linfo('Calling third-party procedure "info", with the same name than another one from basicmod.')`  
!! `end program`  
!!
!! @warning A program using `basicmod` cannot use another module with identical name than those from `basicmod`. The following 
!! piece of code is __incorrect__:  
!! `use basicmod, only: maxpath`  
!! `use module_report !INCORRECT if module_report is a different module from that contained in basicmod` 
!! In particular, this restricction prohibits the use of old `basicmod` versions along with the current one. There are several 
!!possible solutions [here](https://fedoraproject.org/wiki/PackagingDrafts/FortranLibraries#Libraries_with_Interface-Only_Modules), 
!![here](http://stackoverflow.com/questions/10834402/how-to-provide-an-explicit-interface-to-a-library-of-fortran-95-modules-with-i)
!!  and [here](https://software.intel.com/es-es/forums/intel-visual-fortran-compiler-for-windows/topic/520421).  
!! Apart from the obvious solution of renaming old basicmod modules before using them, `basicmod` could be transformed in:  
!! &nbsp;1) __a library of external procedures__, containing modules just for the explicict interfaces; this solution was discarded
!! because it prevents from the use of class atributes.  
!! &nbsp;2) __a library of submodules__; this option will be studied for a future release.  

use module_alloc_bmod
use module_args_bmod
use module_compiler_dependant_bmod
use module_convers_bmod
use module_feed_bmod
use module_files_bmod
use module_math_bmod
use module_os_dependant_bmod
use module_readPVD_bmod
use module_readVTU_bmod
use module_report_bmod
use module_set_bmod
use module_system_bmod
use module_writePVD_bmod
use module_writeVTU_bmod
use module_xml_parser_bmod
use module_xread_bmod
implicit none

end module
