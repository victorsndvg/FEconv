program feconv
!-----------------------------------------------------------------------
! Utility to convert between several mesh and FE field formats
!
! Licensing: This code is distributed under the GNU GPL license.
! Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
! Last update: see variable 'last_update'
!-----------------------------------------------------------------------
use module_convers, only: string
use module_report, only: error
use module_args, only: is_arg
use module_feconv, only: convert
implicit none
character(10) :: last_update = '23/01/2014'

!read and store arguments
if (command_argument_count() == 0) call error('(feconv) command line arguments not found; to show help information: feconv -h')
if (is_arg('-v')) then
  !show version
  print '(a)', 'Utility to convert between several mesh and FE field formats'
  print '(a)', 'Licensing: This code is distributed under the GNU GPL license'
  print '(a)', 'Author: Ingenieria Matematica, http://www.usc.es/ingmat/?lang=en'
  print '(a/)', 'Last update: '//last_update
  stop
elseif (is_arg('-h')) then
  !show help
  print '(a)', 'NAME'
  print '(a)', '    feconv - utility to convert between several mesh and FE field formats'
  print '(a)', ' '
  print '(a)', 'SYNOPSIS'
  print '(a)', '    creamake [OPTIONS] <infile> <outfile>'
  print '(a)', ' '
  print '(a)', 'DESCRIPTION'
  print '(a)', '    The feconv converts the file <infile> into the file <outfile>. The format of each '
  print '(a)', '    file indicated by its extension. Visit the EXAMPLES section of this help to see which formats '
  print '(a)', '    and transformations are available.'
  print '(a)', ' '
  print '(a)', '    The options are as follows:'
  print '(a)', ' '
  print '(a)', '    -h'
  print '(a)', '        Show this help and exit'
  print '(a)', ' '
  print '(a)', '    -v'
  print '(a)', '        Show version information and exit'
  print '(a)', ' '
  print '(a)', '    -is'
  print '(a)', '        Input file contains P2 isoparametrical elements that must be read using Salome ordering (vertices'
  print '(a)', '        and midpoints sandwiched)'
  print '(a)', ' '
  print '(a)', '    -l1'
  print '(a)', '        Output file will contain Lagrange P1 finite elements'
  print '(a)', ' '
  print '(a)', '    -l2'
  print '(a)', '        Output file will contain Lagrange P2 finite elements'
  print '(a)', ' '
  print '(a)', '    -rt'
  print '(a)', '        Output file will contain Raviart-Thomas (face) finite elements of the lowest degree'
  print '(a)', ' '
  print '(a)', '    -nd'
  print '(a)', '        Output file will contain Whitney (edge) finite elements of the lowest degree'
  print '(a)', ' '
  print '(a)', '    -j <opt>'
  print '(a)', '        Ensure positive jacobian in groups with maximal topological dimension. '
  print '(a)', '        <opt> must appear with any of both values, "yes" or "no", in order to force or no the '
  print '(a)', '        positive jacobian. It is recomended to use option "no" with elements of topological '
  print '(a)', '        dimension 2 noncoplanar to a vertical plane. By default <opt> value is set to "yes".'
  print '(a)', ' '
  print '(a)', '    -t <tolerance>'
  print '(a)', '        Tolerance to search mid-points in Lagrange P2 meshes. By default the value '
  print '(a)', '        of tolerance is epsilon for kind real64 (double precission in fortran 77. '
  print '(a)', '        Do not change it unless the conversion fails.'
  print '(a)', ' '
  print '(a)', '    -cm'
  print '(a)', '        Performs a bandwidth optimization using the CutHill-McKee algorithm'
  print '(a)', ' '
  print '(a)', '    -ff'
  print '(a)', '        When a .msh file appears, format applied is FreeFem++ MSH format; if absent, format applied is '
  print '(a)', '        ANSYS MSH format.'
  print '(a)', ' '
  print '(a)', '    -gm'
  print '(a)', '        When a .msh file appears, format applied is Gmsh MSH ASCII format; if absent, format applied is '
  print '(a)', '        ANSYS MSH format.'
  print '(a)', ' '
  print '(a)', '    -es <n>, -e [<n1>,<n2>,<n3>,...]'
  print '(a)', '        Extract subdomain reference <n> as another mesh. When several domain references must be extrated, '
  print '(a)', '        enclose them in square brakets and separate them with commas (without blank spaces!). The name of the '
  print '(a)', '        submesh file is <outfile>. This option is available only for Lagrange P1 meshes.'
  print '(a)', ' '
  print '(a)', '    -p <n>'
  print '(a)', '        When the input mesh is divided in several pieces (like in VTU or MPHTXT formats), <n> is '
  print '(a)', '        the number of the piece to be saved. If not present, when transforming to MFM, all pieces are merged and '
  print '(a)', '        coincident nodes from different pieces are not indentified; use option -glue identify coincident nodes '
  print '(a)', '        from different pieces (not yet implemented)'
  print '(a)', ' '
!  print '(a)', '    -glue'
!  print '(a)', '        When transforming to MFM several pieces, identify coincident nodes from different pieces'
!  print '(a)', '        (not yet implemented)'
!  print '(a)', ' '
!  print '(a)', '    -s <file> (not implemented)'
!  print '(a)', '        A series of fields must be processed; <file> is the name of an ASCII file containing '
!  print '(a)', '        two columns separated by blank spaces: for each line, the first column is a real number'
!  print '(a)', '        representing a parameter value (time, frequency, etc.); the second column is the name'
!  print '(a)', '        of the field file associated with the previous parameter'
!  print '(a)', ' '
  print '(a)', 'FORMATS'
  print '(a)', '    The current available input formats are: ANSYS (.msh), I-Deas Universal (.unv), MD Nastran input file (.bdf),'
  print '(a)', '    Modulef Formatted Mesh (.mfm), Modulef Unformatted Mesh (.mum), FreFem++ Tetrahedral and/or Triangular '
  print '(a)', '    Lagrange P1 Mesh (.msh), FreFem++ Tetrahedral Lagrange P1 Mesh (.mesh), Gmsh MSH ASCII (.msh)'
  print '(a)', ' '
  print '(a)', '    The current available output formats are: VTK-XML Unstructured Grid (.vtu), Modulef Formatted Mesh (.mfm),'
  print '(a)', '    Modulef Unformatted Mesh (.mum), FreFem++ Tetrahedral and/or Triangular Lagrange P1 Mesh (.msh), '
  print '(a)', '    FreFem++ Tetrahedral Lagrange P1 Mesh (.mesh)'
  print '(a)', ' '
  print '(a)', '    For more details, please visit doc/index.xhtml'
  print '(a)', ' '
  print '(a)', 'EXAMPLES'
  print '(a)', '    To transform an ANSYS file into a VTU file:'
  print '(a)', ' '
  print '(a)', '       feconv file1.msh file2.vtu'
  print '(a)', ' '
  print '(a)', '    To transform a MD Nastran input file into Whitney (edge) MFM file preforming Cuthill-McKee optimization:'
  print '(a)', ' '
  print '(a)', '       feconv -nd -cm file1.bdf file2.mfm'
  print '(a)', ' '
  print '(a)', 'AUTHORS'
  print '(a)', '    Written by the members of the Ingenieria Matematica research group (http://www.usc.es/ingmat/?lang=en):'
  print '(a)', '      Iban Constenla,'
  print '(a)', '      Francisco Pena,'
  print '(a)', '      Victor Sande.'
  print '(a)', '    Code related to CutHill-McKee algoritm has been adapted from the one distributed by John Burkardt'
  print '(a)', '    (http://people.sc.fsu.edu/~jburkardt/) under the GNU LGPL license.'
  print '(a)', ' '
  print '(a)', 'REPORTING BUGS'
  print '(a)', '    Report bugs to fran(dot)pena(at)usc(dot)es'
  print '(a)', ' '
  print '(a)', 'DOWNLOAD PAGE'
  print '(a)', '    <http://sourceforge.net/projects/feconv/>'
  print '(a)', ' '
  print '(a)', 'COPYRIGHT'
  print '(a)', '    Copyright © 2010 Free Software Foundation, Inc.  License GPLv3+: GNU GPL version 3 or later'//&
&' <http://gnu.org/licenses/gpl.html>.'
  print '(a)', '    This is free software: you are free to change and redistribute it.  There is NO WARRANTY, '//&
&'to the extent permitted by law.'
  print '(a)', ' '
  stop
end if

!choose converter according to arguments
call convert()

end program
