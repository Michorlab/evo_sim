
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
BUILDDIR = build
OBJS = $(BUILDDIR)/main.o $(BUILDDIR)/MutationHandler.o $(BUILDDIR)/CList.o $(BUILDDIR)/Clone.o $(BUILDDIR)/OutputWriter.o

$(BUILDDIR)/evo_sim : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(BUILDDIR)/evo_sim

$(BUILDDIR)/main.o : main.h main.cpp Clone.h CList.h OutputWriter.h MutationHandler.h
	$(CC) $(CFLAGS) main.cpp -o $(BUILDDIR)/main.o

$(BUILDDIR)/Clone.o : Clone.h Clone.cpp CList.h main.h
	$(CC) $(CFLAGS) Clone.cpp -o $(BUILDDIR)/Clone.o

$(BUILDDIR)/CList.o : CList.h CList.cpp Clone.h main.h
	$(CC) $(CFLAGS) CList.cpp -o $(BUILDDIR)/CList.o

$(BUILDDIR)/OutputWriter.o : OutputWriter.h OutputWriter.cpp main.h Clone.h CList.h
	$(CC) $(CFLAGS) OutputWriter.cpp -o $(BUILDDIR)/OutputWriter.o

$(BUILDDIR)/MutationHandler.o : MutationHandler.h MutationHandler.cpp CList.h main.h
	$(CC) $(CFLAGS) MutationHandler.cpp -o $(BUILDDIR)/MutationHandler.o

CList.h : main.h Clone.h

clean:
	\rm $(BUILDDIR)/*.o $(BUILDDIR)/evo_sim