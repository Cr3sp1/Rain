CC = g++
CFLAGS = -Wall -O3 --std=c++11 
LDFLAGS = -larmadillo
MYLIBS = random.o DataBlocking.o rain.o body.o rainfun.o


main.exe: main.cpp ${MYLIBS}
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

%.o : %.cpp 
	$(CC) -c $< $(CFLAGS) $(LDFLAGS)

run: 
	./main.exe

clean :
	rm *.o main.exe seed.out