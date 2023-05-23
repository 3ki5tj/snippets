#!/usr/bin/env python



''' encode PDB file into a long string '''



import glob, os, sys



def encode(fnpdb):
  ''' encode a single PDB file '''
  s = [ s.strip() for s in open(fnpdb).readlines() ]
  return r'\n'.join( s )



def encodels(ls):
  s = ""
  for fn in ls:
    fn1 = os.path.split( fn )[1]
    name = os.path.splitext( fn1 )[0]
    code = encode(fn)
    s += '"' + name + '": "' + code + '",\n'
  return s



if __name__ == "__main__":
  pdbls = glob.glob("pdb/*.pdb")
  print encodels(pdbls)
