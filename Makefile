# user-editable configuration
CFLAGS  = -march=native -Os
CFLAGS  = -g -Wall -Wextra -Wno-unused

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

# manpages
%.1 : % ; help2man -N ./$^ > $@
manpages: $(BIN) $(MAN)

# bureaucracy
clean   : ; $(RM) $(BIN) $(OBJ)
distclean: clean ; $(RM) $(MAN)
.PHONY  : default clean distclean
