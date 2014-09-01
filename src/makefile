# Libraries ##############################################################

API_DIR =  $(FLUIDPROP)/api/c
INC     = -I$(API_DIR)
LIB     =  ${API_DIR}/libcFPapi.a

#TH_DIR  =  /home/mpini/SOURCE/TDlibrary
#INC_TH  = -I$(TH_DIR)
#LIBTH   =  ${TH_DIR}/TDlibrary.a

# Define compiler and optimizer's flags
LNK    = icpc 
#CMP    = $(LNK) -c
CMP    = $(LNK) -c -traceback
CFLAGS = -DWITH_OPENGL

# List of objects
OBJFLD = fluid_model.o fluid_model_pig.o fluid_model_pvdw.o fluid_model_ppr.o fluid_model_flp.o
OBJTSP = transport_model.o
OBJ    = main.o $(OBJFLD) $(OBJTSP)
OBJ2   = main2.o $(OBJFLD) $(OBJTSP)

# Tasks ##############################################################

all:  TDlibrary.exe SU2_SRT.exe  

.SUFFIXES: .c .cpp 

.c.o:
	$(CMP) $(OPT) $(CFLAGS) $(INC) $(INC_TH) $<
		
.cpp.o:
	$(CMP) $(OPT) $(CFLAGS) $(INC) $(INC_TH) $<

clean:
	rm *.o TDlibrary.exe SU2_SRT.exe
purge:
	rm *.o 


# Executable #########################################################

# ThermoPropPP.a: $(OBJLIB) 
# 	ar cq $@ $(OBJLIB) 

TDlibrary.exe: $(OBJ) 
	$(LNK) $(OBJ) -o TDlibrary.exe $(LIB) -ldl

SU2_SRT.exe: $(OBJ2) 
	$(LNK) $(OBJ2) -o SU2_SRT.exe $(LIB) -ldl
	