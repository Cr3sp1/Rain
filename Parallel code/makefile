CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11 

main.exe : main.o random.o 
	$(CC) random.o main.o -o main.exe $(CFLAGS)
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
run: 
	mpiexec -np 6 ./main.exe 2> error.log
clean :
	rm *.o main.exe 
cleanbtr :
	rm *.btr
cleandat :
	rm *.dat
