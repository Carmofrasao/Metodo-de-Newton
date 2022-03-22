CFLAGS = -Wall -g 
LFLAGS = -lmatheval

objs = main.o utils.o SistLinear.o Metodo_de_Newton_Inexato.o Metodo_de_Newton_Modificado.o Metodo_de_Newton_Padrao.o

all: newtonPC

newtonPC: $(objs)
	gcc $(objs) -o newtonPC $(CFLAGS) $(LFLAGS)
        
main.o: main.c
	cc -c main.c $(CFLAGS)
utils.o: utils.c
	cc -c utils.c $(CFLAGS) 
SistLinear.o: SistLinear.c
	cc -c SistLinear.c $(CFLAGS)  
Metodo_de_Newton_Inexato.o: Metodo_de_Newton_Inexato.c
	cc -c Metodo_de_Newton_Inexato.c $(CFLAGS) 
Metodo_de_Newton_Modificado.o: Metodo_de_Newton_Modificado.c
	cc -c Metodo_de_Newton_Modificado.c $(CFLAGS)
Metodo_de_Newton_Padrao.o: Metodo_de_Newton_Padrao.c
	cc -c Metodo_de_Newton_Padrao.c $(CFLAGS) 
                        
clean:
	-rm -f *.o *~ 
        
purge: clean
	-rm -f newtonPC