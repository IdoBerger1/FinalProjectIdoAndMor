

CC = g++
OBJDIR = Link
EXECUTABLE = OutputProg


OS := $(shell uname)

	MKLROOT = /opt/intel/mkl
	CPPFLAGS = -Ofast -std=c++17 -fopenmp -DMKL_ILP64 -m64 -I${MKLROOT}/include -mavx -mfma -march=native -Wall
 	LINKFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -fopenmp 




all: Release

Release: $(OBJDIR)/main.o $(OBJDIR)/t_timer.o
	$(CC) $(OBJDIR)/main.o $(OBJDIR)/t_timer.o $(LINKFLAGS) -o $(EXECUTABLE)


$(OBJDIR)/main.o: 	main.cpp
	$(CC) $(CPPFLAGS) -c main.cpp -o $(OBJDIR)/main.o 
$(OBJDIR)/t_timer.o: 	t_timer.cpp
	$(CC) $(CPPFLAGS) -c t_timer.cpp -o $(OBJDIR)/t_timer.o


