DEST_LIBS= phaselib.so levelzerolib.so

all: ${DEST_LIBS}

clean:
	rm -rf *.o *.so *.pyc *.py~

install: ${DEST_LIBS}
	cp $^ ${HOME}/lib/

phaselib.so: phaselib.c
	gcc -shared -I${HOME}/x86_64/include -L${HOME}/x86_64/lib -L${NADC_EXTERN}/lib -lSCIA -lNADC -lhdf5 -lhdf5_hl -Wl,-soname,phaselib -o $@ -fPIC $<

levelzerolib.so: levelzerolib.c
	gcc -shared -I${HOME}/x86_64/include -I${NADC_EXTERN}/include -L${HOME}/x86_64/lib -L${NADC_EXTERN}/lib -lSCIA -lNADC -lhdf5 -lhdf5_hl -Wl,-soname,phaselib -o $@ -fPIC $<
