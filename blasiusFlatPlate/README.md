# blasiusFlatPlate — strato limite laminare su lastra piana (Blasius)

Caso OpenFOAM 2D, stazionario, incomprimibile e laminare (`simpleFoam`), pensato
come test rapido dell'installazione e di validazione contro la soluzione di Blasius.

| Parametro | Valore |
|---|---|
| Velocità asintotica U | 1 m/s |
| Viscosità cinematica ν | 1e-5 m²/s |
| Lunghezza lastra L | 1 m → Re_L = 1e5 |
| Dominio | x ∈ [-0.2, 1] m, y ∈ [0, 0.1] m |
| Mesh | 300 × 100 celle (30 000), addensata su bordo d'attacco e parete |

Contorni: `inlet` velocità imposta; `outlet` e `top` pressione fissa (p = 0);
`symmetryFront` (x < 0) simmetria; `plate` (0 ≤ x ≤ 1) parete no-slip;
`frontAndBack` empty (2D).

## Come lanciarlo

```bash
source /path/to/OpenFOAM/etc/bashrc   # se non già attivo
./Allrun                              # blockMesh + checkMesh + simpleFoam + confronto
./Allclean                            # riporta il caso allo stato iniziale
```

`Allrun` usa `simpleFoam` se presente, altrimenti `foamRun -solver incompressibleFluid`
(openfoam.org v11+). Converge in ~450 iterazioni (~20 s su un core).

Su Windows: va eseguito in WSL / Linux (o in un ambiente OpenFOAM per Windows) dopo aver
convertito i fine riga se necessario (`dos2unix Allrun Allclean`).

## Validazione

`plotBlasius.py` (python3 + numpy, matplotlib opzionale) legge mesh e campo `U` dell'ultimo
tempo, confronta i profili a x = 0.25 / 0.5 / 0.75 m con Blasius in variabile di similitudine
η = y·√(U/νx) e salva `blasius_comparison.png`. Risultato atteso:

```
 x [m]     Re_x  max|du/U|    Cf CFD  Cf Blas.  err Cf   d99 CFD d99 Blas.
 0.251    25088     0.0022   0.00416   0.00419   -0.7%    7.61mm    7.78mm
 0.501    50071     0.0014   0.00296   0.00297   -0.4%   10.85mm   10.99mm
 0.751    75054     0.0012   0.00242   0.00242   -0.2%   13.31mm   13.45mm
```

Per ParaView: aprire `blasiusFlatPlate.foam` (creato da `Allrun`). Il campo
`wallShearStress` sulla lastra viene scritto dal function object in `controlDict`.

## Parallelo (opzionale, per la workstation)

```bash
blockMesh
# numberOfSubdomains e simpleCoeffs.n in system/decomposeParDict (default 4 = 4x1x1)
decomposePar
mpirun -np 4 simpleFoam -parallel > log.solver
reconstructPar -latestTime
python3 plotBlasius.py
```
