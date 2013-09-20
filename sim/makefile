#CC=i586-mingw32msvc-g++
#CC=x86_64-w64-mingw32-g++
#CFLAGS=-Wall -Wno-unknown-pragmas -DARCGIS -static-libgcc -static-libstdc++ #-fopenmp

CC=g++ -std=c++11
CFLAGS=-Wall -fopenmp -DNDEBUG #-DARCGIS #-lX11 -pthread #-ltbb

ODIR=obj
PRE_FLAGS=-O3
#PRE_FLAGS=-g

_OBJ = salamander.o main.o utility.o mtbin.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp
	$(CC) $(PRE_FLAGS) -c -o $@ $< $(CFLAGS)

salamander: $(OBJ) obj/main.o
	$(CC) $(PRE_FLAGS) -o salamander.exe $^ $(CFLAGS)
	du -hs ./salamander.exe

clean:
	rm -f $(ODIR)/*.o *~ core salamander.exe
