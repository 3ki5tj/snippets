CC = gcc -O3 -Wall -Wextra
LM = -lm

CX = icc -Wall -Wremarks
LMX =

srcs = $(wildcard test*.c)
bins = $(patsubst %.c,%, $(srcs))
bins_d = $(patsubst %.c,%_d, $(srcs))

deps = $(wildcard *.h) Makefile

all: $(bins_d) $(bins)

$(bins_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bins) : % : %.c $(deps)
	$(CX) -o $@ $< $(LMX)

clean:
	rm -rf $(prog) a.out *~ .*.un~ */*~ */.*.un~ *.dat $(bins) $(bins_d)
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" * ~/Dropbox/code/linalge/
