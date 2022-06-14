# user-editable configuration
CFLAGS  = -march=native -Os
CFLAGS  = -g -Wall -Wextra

# required libraries
LDLIBS  = -lfftw3 -lpng -ljpeg -ltiff

# files
BIN     = krt #krt3d
OBJ     = iio.o

# default target: build all the binaries
default : $(BIN)

# each binary depends on all the object files
$(BIN)  : $(OBJ)

# bureaucracy
clean   : ; $(RM) $(BIN) $(OBJ)
.PHONY  : default clean
