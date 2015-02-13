exe = hdf5test

#Directories for source files, object files and modules
obj = OBJECTS/
src = SRC/
mod = MODULES/

#Fortran compiler
#cff = ifort -O3 -fpe0
cff = ifort -O3
#cff = ifort  -check bounds -fp-stack-check -debug-parameters all -traceback -ftrapuv -g -fpe0
#cff = ifort -g -mcmodel medium -shared-intel -debug -CB -check bounds -check pointer -check uninit -traceback
ff_flag_obj = -c -module $(mod) -I$(mod) -I/opt/local/EOS/nadc_extern/intel64/include 

#C compiler
cc  = icc
#-check bounds -check uninit -debug inline_debug_info -CB -pg
c_flag_obj = -c 
c_flag_lnk = -L$(IFC_LIB) -L$(IFC_MKL_LIB) -lifport -lifcore -lstdc++ -lsvml -L/opt/local/EOS/nadc_extern/intel64/lib -L/opt/local/EOS/nadc_extern/intel64/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -L/usr/lib -lpthread
ff_flag_lnk = -L$(IFC_LIB) -L$(IFC_MKL_LIB) -lifport -lifcore -lstdc++ -lsvml -L/opt/local/EOS/nadc_extern/intel64/lib -L/opt/local/EOS/nadc_extern/intel64/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -L/usr/lib -lpthread

FFOBJS =\
	$(obj)mask_reader.o\
	$(obj)noise_reader.o\
	$(obj)smr_reader.o\
	$(obj)main.o

all: $(exe)

# if you want to mix c(++) and fortran
#COBJS = $(obj)main.o
#$(exe): $(obj) $(mod) $(FFOBJS) $(COBJS)
#	$(cc) $(FFOBJS) $(COBJS) $(c_flag_lnk) $(FFLIBS) -o $(exe)

# if you have everything fortran (well.. you can still link intel c++ compiled stuff)
$(exe): $(obj) $(mod) $(FFOBJS)
	$(cff) $(FFOBJS) $(ff_flag_lnk) $(FFLIBS) -o $(exe)

$(obj):
	mkdir $(obj)    

$(mod):
	mkdir $(mod)    

$(obj)main.o: $(src)main.f90
	$(cff) $(ff_flag_obj) $(src)main.f90 -o $(obj)main.o

$(obj)noise_reader.o:  $(src)noise_reader.f90
	$(cff) $(ff_flag_obj) $(src)noise_reader.f90 -o $(obj)noise_reader.o
     
$(obj)mask_reader.o:  $(src)mask_reader.f90
	$(cff) $(ff_flag_obj) $(src)mask_reader.f90 -o $(obj)mask_reader.o
     
$(obj)smr_reader.o:  $(src)smr_reader.f90
	$(cff) $(ff_flag_obj) $(src)smr_reader.f90 -o $(obj)smr_reader.o
     
%.o: %.f
	$(cff) $(FFLAGS) $<

%.o: %.f90
	$(cff) $(FFLAGS) $<

run: $(exe)
	./$(exe)

clean:
	rm -rf $(obj)
	rm -rf $(mod)
	rm -f $(exe)
