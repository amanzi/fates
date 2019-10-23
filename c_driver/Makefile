#HOME := $(shell /bin/pwd)
 
CC          = icc
F77	    = ifort



##################################################################

#CFLAGS      = -O3 -Wno-deprecated
CFLAGS      = -g -Wno-deprecated
FFLAGS      = 

CLINKER     = $(CC) 
LFLAGS      = $(CFLAGS)  


##################################################################
#Directories
SYSHOME  =  $(HOME)

INCLDIR  = 


LIBDIRS  = -L$(ATSFATES_DIR)/install/lib/ -L$(NETCDF_DIR)/lib/ -L$(MPI_HOME)/lib -L${PNETCDF_DIR}/lib/ -L${PIO_DIR}/lib/


##################################################################
#Libraries
LIBS = -lparteh -lfates -lshared  -lparteh -lpiof -lpioc -lpnetcdf -lnetcdf -lnetcdff  -lgptl  -lhdf5_hl -lhdf5 -lz -lmpi_mpifh -lmpi -lifcore -lm



##################################################################
#Programs
EXECS  = Test


##################################################################
.SUFFIXES: .c .o
.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o $(INCLDIR)

.f.o:
	$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLDIR)


OBJS =  c_driver.o

##################################################################
all:   $(EXECS) 


$(EXECS): $(OBJS)
	$(CLINKER) $(LFLAGS) $(LIBDIRS) -o $(EXECS) $(OBJS) $(LIBS)

clean:
	@rm -f $(OBJS) Test.exe *~ 


##################################################################
Test : Makefile ../install/lib/libfates.a ../install/lib/libshared.a


##################################################################

