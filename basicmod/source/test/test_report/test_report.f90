program test_report
use basicmod, only: report_option, report, debug, info, warning, error  
implicit none

print*,' '
print*, 'call report_option(''file'', ''out.txt'')'
         call report_option( 'file',   'out.txt')
print*, 'call report_option(''folder'', ''out'')'
         call report_option( 'folder',   'out')

print*, 'call report_option(''debug'', ''none'')'
         call report_option( 'debug',   'none')
print*, 'call report(''debug'', ''report debug'')'
         call report( 'debug',   'report debug')
print*, 'call debug(''sub debug'')'
         call debug( 'sub debug')

print*, 'call report_option(''info'', ''std'')'
         call report_option( 'info',   'std')
print*, 'call report(''info'', ''report info'')'
         call report( 'info',   'report info')
print*, 'call info(''sub info'')'
         call info( 'sub info')

print*, 'call report_option(''warning'', ''file'')'
         call report_option( 'warning',   'file')
print*, 'call report(''warning'', ''report warning'')'
         call report( 'warning',   'report warning')
print*, 'call warning(''sub warning'')'
         call warning( 'sub warning')

print*, 'call report_option(''error'', ''all'')'
         call report_option( 'error',   'all')
!print*, 'call report(''error'', ''report error'')'
!         call report( 'error',   'report error')
print*, 'call error(''sub error'')'
         call error( 'sub error')
end program
