subdirs = pca kmeans lj wham

clean:
	for dir in $(subdirs) ; do $(MAKE) -C $$dir clean ; done
	rm -rf $(prog) a.out *~ .*.un~ */*~ */.*.un~ */*/*~ */*/.*.un~
	rstrip.py -Rlv

exclude = --exclude=".*" --exclude="*~"

Dropbox: clean
	rsync -avzL $(exclude) * ~/Dropbox/code/
