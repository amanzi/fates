#HOME := $(shell /bin/pwd)
 
CC          = mpicc
F77	    = mpif90



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


LIBDIRS  = -L../install/lib/ -L$(NETCDF_DIR)/lib/ -L$(MPI_HOME)/lib -L$(HOME)/Coastal/ParallelIO/install/lib/ -L${HOME}/Coastal/pnetcdf-1.11.2/install/lib


##################################################################
#Libraries
LIBS = -lfates -lshared  -lpiof -lpioc -lpnetcdf -lnetcdf -lnetcdff -lgptl -lmpi_mpifh -lmpi -lifcore -lm



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

