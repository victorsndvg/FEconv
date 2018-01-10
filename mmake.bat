creamake configure/creamake.%1.xml
move     Makefile  Makefile.%1
mingw32-make -f    Makefile.%1 clean
mingw32-make -f    Makefile.%1
