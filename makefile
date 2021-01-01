CC       = gcc 
CFLAGS   = -O3 -Wall 
LIBS      = -lXi -lXmu -lglut -lGLEW -lGLU -lm -lGL
OBJDIR   = ../linalglib
OBJS     = $(OBJDIR)/linalglib.o maze.o initShader.o

mazesolver: mazesolver.c $(OBJS)
	$(CC) -o mazesolver mazesolver.c $(OBJS) $(CFLAGS) $(LIBS)

testmaze: testmaze.c $(OBJS)
	$(CC) -o testmaze testmaze.c $(OBJS) $(CFLAGS) $(LIBS)

$(OBJDIR)/%.o: %.c %.h
	$(CC) -c @< -o $@ $(CFLAGS)

