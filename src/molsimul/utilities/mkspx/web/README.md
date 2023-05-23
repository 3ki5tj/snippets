# mkspx web app

psfgen:

```psfgen
package require psfgen
topology top_all27_prot_lipid.inp
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD
segment U {pdb out0.pdb}
coordpdb out0.pdb U
guesscoord
writepdb out.pdb
writepsf out.psf
```
