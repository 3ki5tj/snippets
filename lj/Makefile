CC = gcc -O3 -Wall -Wextra
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
bin2ds = $(patsubst %,%2d, $(bins))
deps = $(wildcard *.h)

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

$(bin2ds): %2d : %.c $(deps)
	$(CC) -DD=2 $< -o $@ $(LM)

clean:
	rm -rf $(bins) $(bin2ds) a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" * ~/Dropbox/code/kmeans/
