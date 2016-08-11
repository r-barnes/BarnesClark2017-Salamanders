#For 'march' options see: https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/i386-and-x86-64-Options.html
#Itasca: -march=nehalem
CC=g++
CFLAGS=-Wall --std=c++11 -ffast-math -march=native -fopenmp  #-DNDEBUG

ODIR=obj
PRE_FLAGS=-O3 -g

_OBJ = salamander.o mtbin.o temp.o phylo.o random.o simulation.o params.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp
	$(CC) $(PRE_FLAGS) -c -o $@ $< $(CFLAGS)

salamander: $(OBJ) obj/main.o
	$(CC) $(PRE_FLAGS) -o salamander.exe $^ $(CFLAGS)
	du -hs ./salamander.exe

test: $(OBJ) obj/test.o
	$(CC) $(PRE_FLAGS) -o test.exe $^ $(CFLAGS)
	du -hs ./test.exe	

clean:
	rm -f $(ODIR)/*.o *~ core salamander.exe
