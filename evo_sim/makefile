CC = g++ -std=c++11 -static-libstdc++ -lpthread
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
BUILDDIR = build
OBJS = $(BUILDDIR)/main.o $(BUILDDIR)/MutationHandler.o $(BUILDDIR)/CList.o $(BUILDDIR)/Clone.o $(BUILDDIR)/OutputWriter.o

$(shell   mkdir -p $(BUILDDIR))

$(BUILDDIR)/evo_sim : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(BUILDDIR)/evo_sim

$(BUILDDIR)/main.o : main.cpp Clone.h CList.h OutputWriter.h MutationHandler.h main.h 
	$(CC) $(CFLAGS) main.cpp -o $(BUILDDIR)/main.o

$(BUILDDIR)/Clone.o : Clone.cpp Clone.h CList.h OutputWriter.h MutationHandler.h main.h
	$(CC) $(CFLAGS) Clone.cpp -o $(BUILDDIR)/Clone.o

$(BUILDDIR)/CList.o : CList.cpp Clone.h CList.h OutputWriter.h MutationHandler.h main.h
	$(CC) $(CFLAGS) CList.cpp -o $(BUILDDIR)/CList.o

$(BUILDDIR)/OutputWriter.o : OutputWriter.cpp Clone.h CList.h Clone.h CList.h OutputWriter.h MutationHandler.h main.h
	$(CC) $(CFLAGS) OutputWriter.cpp -o $(BUILDDIR)/OutputWriter.o

$(BUILDDIR)/MutationHandler.o : MutationHandler.cpp Clone.h CList.h OutputWriter.h MutationHandler.h main.h
	$(CC) $(CFLAGS) MutationHandler.cpp -o $(BUILDDIR)/MutationHandler.o

CList.h : main.h Clone.h

clean:
	\rm $(BUILDDIR)/*.o $(BUILDDIR)/evo_sim