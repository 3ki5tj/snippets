Comparison of RMSD

```
make polymer && ./polymer -t100000000 --Kd=0.1 --his=rad_Langevin.his
make polymer && ./polymer -t100000000 --Kd=0.1 --th=vrescale --his=rad_vrescale.his
```

```
plot "rad_Langevin.his" u 1:2 w l, "rad_vrescale.his" u 1:($2*$1**3/26.5) w l
```
