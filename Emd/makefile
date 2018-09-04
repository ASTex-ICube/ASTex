CC = gcc
CFLAGS = -g -lm
OBJECTS = main.o emd.o

main : $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) -o main

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o main
