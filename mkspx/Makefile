clean:
	rm -f *~ .*~ */*~ */.*~ out* topol.top posre.itp \#*
	rstrip.py * js/*.js

bdir=/Bossman/cz1/snippets/mkspx

Bossman:
	rm -f *~ js/*~ .*.un~ js/.*.un~
	rstrip.py -R
	../htmlpack.py
	rsync -avzL js/*.js $(bdir)/js
	rsync -avz mkspx*.html *.inp $(bdir)/
