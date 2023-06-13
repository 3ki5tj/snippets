#CC = icc -g -O2 -Wall -Wremarks
#LM =
CC = gcc -g -O3 -Wall -Wextra
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard ../*.h) $(wildcard *.h) Makefile
includes = -I../../../../common/c \
           -I../../../systems/lj/c \
           -I../../../basic/md/c \
           -I../../../methods

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $(includes) $< -o $@ $(LM)

clean:
	rm -rf $(bins) a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos

