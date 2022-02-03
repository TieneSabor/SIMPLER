PROJNAME := Final
CC := gcc
SRC := $(wildcard src/*.c)
NOD := $(notdir $(SRC)) #SRC w/o dir.
OBJ := $(patsubst %.c,bin/%.o,$(NOD))
INC := -I inc/
CFLAGS := -g -Wa,--warn -ffunction-sections -fdata-sections
LINKFLAGS := -Wl,--gc-sections
MATHFLAGS := -lm
PREPRO := -D BASIC_LINEAR

prt:
	@echo $(OBJ)

all: $(OBJ)
	$(CC) $(PREPRO) -o bin/$(PROJNAME) $(OBJ) $(LINKFLAGS) $(MATHFLAGS)

bin/%.o: src/%.c
	$(CC) $(PREPRO) $(CFLAGS) $(INC) -c $< -o $@ $(MATHFLAGS)

clean:
	rm bin/*
