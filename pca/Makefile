CC = gcc -O3 -Wall -Wextra
LM = -lm

prog = ho

ho: ho.c mtrand.h eig.h lu.h
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(prog) a.out *~ .*.un~ */*~ */.*.un~
	rstrip.py
