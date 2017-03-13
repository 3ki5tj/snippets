clean:
	rm -f *~ .*~ */*~ */.*~ out* topol.top posre.itp \#*
	rstrip.py * js/*.js

bdir=/Bossman/cz1/snippets/ewald

Bossman:
	rm -f *~ js/*~ .*.un~ js/.*.un~
	rstrip.py -R js/* *
	../htmlpack.py
	rsync -avzL js/*.js $(bdir)/js
	rsync -avz ion.html $(bdir)/
