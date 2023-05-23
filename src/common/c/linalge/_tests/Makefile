INCLUDES = -I..
CC = gcc -O3 -g -Wall -Wextra $(INCLUDES)
LM = -lm

CC_DBG = gcc -O0 -g -pg -Wall -Wextra $(INCLUDES)
LM_DBG = -lm

srcs = $(wildcard test*.c)
bins = $(patsubst %.c,%, $(srcs))
bins_d = $(patsubst %.c,%_d, $(srcs))

deps = $(wildcard *.h) Makefile

all: $(bins)

$(bins_d) : %_d : %.c $(deps)
	$(CC_DBG) -o $@ $< $(LM_DBG)

$(bins) : % : %.c $(deps)
	$(CC) -o $@ $< $(LM)

clean:
	rm -rf $(prog) a.out *~ .*.un~ */*~ */.*.un~ *.dat $(bins) $(bins_d)
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" * ~/Dropbox/code/linalge/
