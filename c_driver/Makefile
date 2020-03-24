#HOME := $(shell /bin/pwd)
 

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


LIBDIRS  = -L$(TPLS_DIR)/lib/ -L$(TPLS_DIR)/lib64/ -L${MPIHOME}/lib/


#################################################################q#Libraries
LIBS = -lfates -lshared  -lparteh -lpiof -lpioc -lpnetcdf -lnetcdf -lnetcdff  -lhdf5_hl -lhdf5 -lz -lmpi_mpifh -lifcore -lm



##################################################################
#Programs
EXECS  = Test


##################################################################
.SUFFIXES: .c .o
.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o $(INCLDIR)

.f.o:
	$(FC) $(FFLAGS) -c $*.f -o $*.o $(INCLDIR)


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

