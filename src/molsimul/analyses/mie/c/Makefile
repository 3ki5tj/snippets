INCLUDES = -I../../../../common/c/util \
           -I../../../../common/c/rand \
           -I../../../../common/c/ave \
           -I../../../../common/c/argopt

CC = gcc -O3 -Wall -Wextra $(INCLUDES)
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard *.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos vgcore.*
	rstrip.py -Rlv

zip:
	zip mie.zip *.c *.h Makefile \
	    doc/Makefile doc/*.tex doc/*.nb
