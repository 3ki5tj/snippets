CC = icc -g -O2 -Wall -Wremarks
LM =
#CC = gcc -g -O3 -Wall -Wextra
#LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard ../*.h) $(wildcard *.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" * ~/Dropbox/code/wham/
