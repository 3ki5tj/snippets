INCLUDES = -I../../../common/c/vct \
           -I../../../common/c/rand \
           -I../../systems/ising2d

CC = gcc -O3 -Wall -Wextra $(INCLUDES)
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard *.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" * ~/Dropbox/code/ising/
