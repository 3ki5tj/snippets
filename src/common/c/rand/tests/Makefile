CC = gcc -O3 -Wall -Wextra
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard *.h)
includes = -I.. -I../..

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $(includes) $< -o $@ $(LM)

clean:
	rm -rf $(bins) a.out *~ .*.un~ */*~ */.*.un~ *.dat
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" * ~/Dropbox/code/mtrand/

