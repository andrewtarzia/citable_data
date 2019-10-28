In this folder are the original input files for Gaussian16 and the final output of the frequency calculation as .log (Gaussian output from which values were extracted), .pdb/.xyz (final optimised structures).

Some notes:
-----------

* For 4_D_dmso:
    * The convergence with the default solvent mesh was poor, so I used surface=SAS and optimised the structure (this is not shown in the original .gjf). 
    * Once optimised, I then took the output and optimised that using the same setup as all other systems.

* Long optimisations were restarted using:
    * geom=(checkpoint,step=X) if energy of step X < energy of the latest step.
    * calcfc (if following an optimisation) and rcfc (if following a test ferquency calculation) to calc and/or read force constants at the start of the optimisation.

* All frequency calculations were run with 'freq'.

* default method line:
    * #P B3LYP/GenECP opt=(Tight,steep) SCF=(YQC,MaxCycle=900) int=(Grid=Ultrafine) EmpiricalDispersion=GD3 SCRF=(PCM,Solvent=Acetonitrile)

* Solvents:
    * mecn = Acetonitrile
    * dmso = DiMethylSulfoxide

* default ECP section:

```
C H N 0
6-31G(d)
****
Pd 0
SDD
****

Pd 0
SDD
```
