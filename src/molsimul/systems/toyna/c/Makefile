INCLUDES = -I../../../../common/c/vct \
           -I../../../../common/c/rand \
           -I../../../basic/md/c

CC = gcc -O3 -Wall -Wextra $(INCLUDES)
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard *.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) $(bin2ds) a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" --exclude=".*" \
	  * ~/Dropbox/code/toyna/
