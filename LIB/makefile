# config file.
include ../config/config.inc

src =  isnumb.f lowor.f

obj = $(src:.f=.o)

libsmarddabit$(DBG).a: $(obj)  randlx.o
	ar ru libsmarddabit$(DBG).a $(obj) randlx.o
	ranlib libsmarddabit$(DBG).a

randlx.o:
	gfortran -g -fbacktrace -fbounds-check -fexceptions -c randlx.f

tar: $(src) README makefile
	tar cvf fort.tar $(src)  README makefile

clean: 
	rm -f $(obj) 

libclean: clean
	rm -f libsmarddabit$(DBG).a
