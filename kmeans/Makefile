CC = gcc -O3 -Wall -Wextra
LM = -lm

bins = test

bins : %: %.c kmeans.h cholesky.h
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(prog) a.out *~ .*.un~ */*~ */.*.un~ *.dat $(bins)
	rstrip.py -Rlv

Dropbox: clean
	rsync -avzL --exclude="*~" * ~/Dropbox/code/kmeans/
