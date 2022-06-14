# user-editable configuration
CFLAGS  = -g -Wall -Wextra
CFLAGS  = -march=native -O3

# required libraries
LDLIBS  = -lfftw3 -lpng -ljpeg -ltiff

# files
BIN     = krt #krt3d
OBJ     = iio.o
MAN     = krt.1

# default target: build all the binaries
default : $(BIN)

# each binary depends on all the object files
$(BIN)  : $(OBJ)

# build manpages from --help output
%.1 : % ; help2man -N ./$^ > $@
manpages: $(BIN) $(MAN)

# bureaucracy
clean   : ; $(RM) $(BIN) $(OBJ)
distclean: clean ; $(RM) $(MAN)
.PHONY  : default clean distclean
