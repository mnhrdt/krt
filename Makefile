# user-editable configuration
CFLAGS  = -g -Wall -Wextra
CFLAGS  = -march=native -O3 -fopenmp

# required libraries
LDLIBS  = -lm -lfftw3 -lpng -ljpeg -ltiff

# files
BIN     = krt #krt3d
OBJ     = iio.o
MAN     = krt.1

# default target: build all the binaries
default : $(BIN)

# each binary depends on all the object files
$(BIN)  : $(OBJ)

# build manpages from --help output
%.1 : % ; SOURCE_DATE_EPOCH=1666666666 help2man -N ./$^ > $@
manpages: $(BIN) $(MAN)

# bureaucracy
clean   : ; $(RM) $(BIN) $(OBJ)
distclean: clean doc-clean ; $(RM) $(MAN)
.PHONY  : default clean distclean manpages doc

# build latex article and associated experiments
doc     : $(BIN) ; $(MAKE) -C article ; $(MAKE) -C slides
doc-clean : ; $(MAKE) -C article clean ; $(MAKE) -C slides clean
