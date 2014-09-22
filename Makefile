all: phaselib.so

clean:
	rm -rf *.o *.so *.pyc *.py~

phaselib.so: phaselib.c
	gcc -shared -I${HOME}/x86_64/include -L${HOME}/x86_64/lib -L${NADC_EXTERN}/lib -lSCIA -lNADC -lhdf5 -lhdf5_hl -Wl,-soname,phaselib -o $@ -fPIC $<
