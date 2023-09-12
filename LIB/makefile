# config file.
include ../config/config.inc

src =  isnumb.f lowor.f

obj = $(src:.f=.o)

libsmarddabit$(DBG).a: $(obj) 
	ar ru libsmarddabit$(DBG).a $(obj) 
	ranlib libsmarddabit$(DBG).a

tar: $(src) README makefile
	tar cvf fort.tar $(src)  README makefile

clean: 
	rm -f $(obj) 

libclean: clean
	rm -f libsmarddabit$(DBG).a
