CC = gcc -O3 -Wall -Wextra
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
bins2d = $(patsubst %.c,%_2d, $(srcs))
deps = $(wildcard *.h) Makefile

all: $(bins) $(bins2d)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

$(bins2d): %_2d : %.c $(deps)
	$(CC) -DD=2 $< -o $@ $(LM)

clean:
	rm -rf $(bins) $(bins2d) a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos *.xyz
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" --exclude=".*" \
	  * ~/Dropbox/code/polymer/
