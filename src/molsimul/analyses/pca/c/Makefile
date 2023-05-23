INCLUDES = -I../../../../common/c/rand \
           -I../../../../common/c/linalge

CC = gcc -O3 -Wall -Wextra $(INCLUDES)
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c, %, $(srcs))
deps = $(wildcard *.h)

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) a.out *~ .*.un~ */*~ */.*.un~
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" * ~/Dropbox/code/pca/
