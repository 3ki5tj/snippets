CC = gcc -O3 -Wall -Wextra
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard *.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) a.out *~ .*.un~ */*~ */.*.un~ *.pos vgcore.* *.his
	rstrip.py -Rlv

zip:
	zip mie.zip *.c *.h Makefile
