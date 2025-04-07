CC=qcc
PREFLAGS=-fopenmp -Wall -O2 -DLAYERS=1 -g -disable-dimensions
POSTFLAGS= -L$(BASILISK)/ppr -lppr -lm
TARGET=river
DEPS=
DIRS=plots logs

$(TARGET) : $(TARGET).o $(DEPS)
	$(CC) $(PREFLAGS) -o $(TARGET) $(TARGET).o $(DEPS) $(POSTFLAGS)

$(TARGET).o : $(TARGET).c
	$(CC) $(PREFLAGS) -c -I. $(TARGET).c $(POSTFLAGS)


$(shell mkdir -p $(DIRS))

clean:
	$(RM) $(TARGET) $(TARGET).o $(DEPS)
	rm -rf plots/ logs/

remake:
	$(RM) $(TARGET) $(TARGET).o $(DEPS)
	make

