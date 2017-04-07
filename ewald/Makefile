clean:
	rm -f *~ .*~ */*~ */.*~ out* topol.top posre.itp \#*
	rstrip.py *.html *.c *.h README js/*.js

bdir=/Bossman/cz1/snippets/ewald

Bossman:
	$(MAKE) clean
	../htmlpack.py
	rsync -avzL js/*.js $(bdir)/js
	rsync -avz ion.html $(bdir)/
