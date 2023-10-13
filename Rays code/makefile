CC = g++
CFLAGS = -Wall -O3 --std=c++11 
MYLIBS = body.o ray.h VectorOperations.h VectorStat.h


main.exe: main.cpp ${MYLIBS}
	$(CC) -o $@ $^ $(CFLAGS) 

%.o : %.cpp 
	$(CC) -c $< $(CFLAGS)

run: 
	./main.exe

clean :
	rm *.o main.exe