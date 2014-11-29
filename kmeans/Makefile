CC = gcc -O3 -Wall -Wextra
LM = -lm

prog = ho

test: test.c kmeans.h cholesky.h
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(prog) a.out *~ .*.un~ */*~ */.*.un~ *.dat
	rstrip.py
