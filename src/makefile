COMPILER = g++
FLAGS = -std=c++14 -pipe -Wall -Wextra -Wshadow -O3 -ffast-math -flto -march=native -pthread -DNDEBUG
OBJS = main.o util.o hssp.o hypervolume.o test.o nth_subsetsum.o
OBJS += hv-2.0rc2-src/hv.o
OBJS += HVC/hvc.o HVC/avl.o HVC/hvc-class.o HVC/io.o

all:
	$(MAKE) -C hv-2.0rc2-src/
	$(MAKE) -C HVC/
	$(MAKE) hvsubset

hvsubset: $(OBJS)
	$(COMPILER) $(FLAGS) -o hssp $(OBJS)

main.o: main.cpp hssp.h test.h
	$(COMPILER) $(FLAGS) -c main.cpp

util.o: util.cpp util.h
	$(COMPILER) $(FLAGS) -c util.cpp

hssp.o: hssp.cpp hssp.h hypervolume.h util.h nth_subsetsum.h HVC/hvc-class.c HVC/hvc-class.h
	$(COMPILER) $(FLAGS) -c hssp.cpp

hypervolume.o: hypervolume.cpp hypervolume.h util.h hv-2.0rc2-src/hv.c hv-2.0rc2-src/hv.h HVC/hvc.c HVC/hvc.h HVC/avl.c HVC/avl.h HVC/hvc-class.c HVC/hvc-class.h HVC/io.c HVC/io.h
	$(COMPILER) $(FLAGS) -c hypervolume.cpp

test.o: test.cpp test.h hssp.h util.h
	$(COMPILER) $(FLAGS) -c test.cpp

nth_subsetsum.o: nth_subsetsum.cpp nth_subsetsum.h
	$(COMPILER) $(FLAGS) -c nth_subsetsum.cpp

clean:
	rm -f hvsubset *.o
	$(MAKE) -C hv-2.0rc2-src/ clean
	$(MAKE) -C HVC/ clean
