# Structural results and literature comparison

Campaign closed; existing energies only. Numerical release remains pending.

## Scope and conventions

Ten cubic solids: AlN, AlP, BN, BP, C, LiCl, MgO, MgS, Si, SiC. LiF and LiH are excluded consistently across methods.
352 unique energies, 36 independent curves and 40 method/solid rows. AlN/AlP use v01-v09; others v01-v10.
BN and C reuse the AE curve in both hybrid representations. The hybrid label denotes mixed AE/GTH treatment, not a hybrid exchange-correlation functional.
All signed errors below use calculation minus the same Goldzak experimental reference, corrected to the static-lattice convention. No method-specific reference correction is fitted.
Results with different source protocols are kept separate. Each comparison uses the exact common solid set, printed below.
This is a comparison with the retrievable relevant benchmark literature, not an exhaustive inventory of every published value.

## Native results

### Lattice constant a0 (angstrom)

| Solid | Experiment | GX | AE | Hybrid direct | Hybrid +1c |
| --- | --- | --- | --- | --- | --- |
| AlN | 4.3680 | 4.2811 | 4.3025 | 4.2747 | 4.2746 |
| AlP | 5.4480 | 5.3622 | 5.3781 | 5.3616 | 5.3621 |
| BN | 3.5930 | 3.5777 | 3.5675 | 3.5675 | 3.5675 |
| BP | 4.5250 | 4.4439 | 4.4696 | 4.4527 | 4.4527 |
| C | 3.5530 | 3.5345 | 3.5292 | 3.5292 | 3.5292 |
| LiCl | 5.0720 | 4.9081 | 4.9155 | 4.9178 | 4.9179 |
| MgO | 4.1890 | 4.1059 | 4.0935 | 4.0923 | 4.0919 |
| MgS | 5.1880 | 5.0689 | 5.0536 | 5.0700 | 5.0697 |
| Si | 5.4210 | 5.3386 | 5.3740 | 5.3383 | 5.3383 |
| SiC | 4.3470 | 4.2517 | 4.2920 | 4.2556 | 4.2556 |

### Bulk modulus B0 (GPa)

| Solid | Experiment | GX | AE | Hybrid direct | Hybrid +1c |
| --- | --- | --- | --- | --- | --- |
| AlN | 206.0000 | 239.0714 | 225.8898 | 240.4166 | 240.4005 |
| AlP | 87.0000 | 94.3363 | 101.3488 | 96.3581 | 93.1207 |
| BN | 388.5000 | 423.3935 | 420.8738 | 420.8738 | 420.8738 |
| BP | 176.5000 | 187.9280 | 183.6395 | 185.8573 | 185.8553 |
| C | 453.3000 | 457.3333 | 491.3550 | 491.3550 | 491.3550 |
| LiCl | 37.3000 | 47.9760 | 46.3654 | 47.1319 | 47.1640 |
| MgO | 173.0000 | 188.2759 | 198.7086 | 213.1266 | 207.3747 |
| MgS | 81.0000 | 87.8957 | 90.2744 | 84.1736 | 86.2197 |
| Si | 100.3000 | 92.7635 | 100.5632 | 94.1294 | 94.1297 |
| SiC | 228.9000 | 263.8872 | 241.6757 | 263.7281 | 263.7096 |

## Goldzak2022

[Primary source](https://doi.org/10.1063/5.0119633)

### a_A: matched errors

Solids: AlN, AlP, BN, BP, C, LiCl, MgO, MgS, Si, SiC.

| Method | N | ME | MAE | RMSE | MaxAE | MARE (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Skala GAPW_XC-GTH | 10 | -0.08313 | 0.08313 | 0.09267 | 0.16391 | 1.75251 |
| Skala GAPW-AE | 10 | -0.07284 | 0.07284 | 0.08389 | 0.15650 | 1.54734 |
| Skala hybrid direct | 10 | -0.08443 | 0.08443 | 0.09212 | 0.15417 | 1.79501 |
| Skala hybrid +1c | 10 | -0.08445 | 0.08445 | 0.09215 | 0.15414 | 1.79564 |
| HF | 10 | 0.05190 | 0.05630 | 0.07880 | 0.18100 | 1.11788 |
| MP2 | 10 | -0.00730 | 0.02010 | 0.02514 | 0.05100 | 0.43588 |
| SCS-MP2 | 10 | 0.00770 | 0.01250 | 0.01563 | 0.03500 | 0.28089 |
| SOS-MP2 | 10 | 0.01490 | 0.01490 | 0.01786 | 0.03400 | 0.32534 |

### a_A: literature values

| Solid | HF | MP2 | SCS-MP2 | SOS-MP2 |
| --- | --- | --- | --- | --- |
| AlN | 4.3650 | 4.3880 | 4.3890 | 4.3890 |
| AlP | 5.5420 | 5.4440 | 5.4650 | 5.4750 |
| BN | 3.5960 | 3.5960 | 3.6030 | 3.6060 |
| BP | 4.5840 | 4.4950 | 4.5170 | 4.5280 |
| C | 3.5470 | 3.5400 | 3.5500 | 3.5540 |
| LiCl | 5.2530 | 5.0210 | 5.0590 | 5.0780 |
| MgO | 4.1760 | 4.2270 | 4.2240 | 4.2230 |
| MgS | 5.2810 | 5.1710 | 5.1910 | 5.2010 |
| Si | 5.5080 | 5.3990 | 5.4250 | 5.4370 |
| SiC | 4.3710 | 4.3500 | 4.3580 | 4.3620 |

### B0_GPa: matched errors

Solids: AlN, AlP, BN, BP, C, LiCl, MgO, MgS, Si, SiC.

| Method | N | ME | MAE | RMSE | MaxAE | MARE (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Skala GAPW_XC-GTH | 10 | 15.10609 | 16.61339 | 20.46400 | 34.98723 | 10.95969 |
| Skala GAPW-AE | 10 | 16.88941 | 16.88941 | 20.34974 | 38.05495 | 10.33794 |
| Skala hybrid direct | 10 | 20.53504 | 21.76915 | 26.12247 | 40.12665 | 12.43324 |
| Skala hybrid +1c | 10 | 19.84030 | 21.07436 | 25.22016 | 38.05495 | 11.98812 |
| HF | 10 | 12.53000 | 15.23000 | 21.79979 | 47.30000 | 7.75101 |
| MP2 | 10 | 1.33000 | 5.71000 | 7.59230 | 15.20000 | 3.20228 |
| SCS-MP2 | 10 | -0.56000 | 4.10000 | 5.61017 | 13.50000 | 2.44803 |
| SOS-MP2 | 10 | -1.46000 | 3.88000 | 4.99520 | 12.60000 | 2.62064 |

### B0_GPa: literature values

| Solid | HF | MP2 | SCS-MP2 | SOS-MP2 |
| --- | --- | --- | --- | --- |
| AlN | 226.1000 | 200.0000 | 201.3000 | 201.9000 |
| AlP | 93.9000 | 93.2000 | 92.0000 | 91.4000 |
| BN | 430.2000 | 396.5000 | 394.0000 | 392.7000 |
| BP | 174.8000 | 180.9000 | 176.2000 | 174.1000 |
| C | 500.6000 | 466.8000 | 460.1000 | 456.8000 |
| LiCl | 29.6000 | 38.2000 | 36.5000 | 35.7000 |
| MgO | 180.2000 | 157.8000 | 159.5000 | 160.4000 |
| MgS | 76.9000 | 83.0000 | 81.4000 | 80.7000 |
| Si | 102.1000 | 99.6000 | 98.2000 | 97.6000 |
| SiC | 242.7000 | 229.1000 | 227.0000 | 225.9000 |

## Zhang2018

[Primary source](https://pure.mpg.de/rest/items/item_2599505_10/component/file_2618826/content)

SI Tables III/V: use the Uncorr. (static electronic) columns, not Corr. columns that include zero-point motion. SCAN was evaluated on PBE orbitals/densities. AlN is absent.

### a_A: matched errors

Solids: AlP, BN, BP, C, LiCl, MgO, MgS, Si, SiC.

| Method | N | ME | MAE | RMSE | MaxAE | MARE (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Skala GAPW_XC-GTH | 9 | -0.08272 | 0.08272 | 0.09330 | 0.16391 | 1.72627 |
| Skala GAPW-AE | 9 | -0.07366 | 0.07366 | 0.08568 | 0.15650 | 1.55254 |
| Skala hybrid direct | 9 | -0.08344 | 0.08344 | 0.09198 | 0.15417 | 1.75701 |
| Skala hybrid +1c | 9 | -0.08346 | 0.08346 | 0.09202 | 0.15414 | 1.75766 |
| LDA | 9 | -0.03244 | 0.03244 | 0.04302 | 0.10600 | 0.68899 |
| PBE | 9 | 0.04511 | 0.04511 | 0.04887 | 0.07900 | 0.96767 |
| PBEsol | 9 | 0.00733 | 0.01133 | 0.01352 | 0.02500 | 0.24954 |
| M06-L | 9 | 0.00478 | 0.01278 | 0.01772 | 0.04500 | 0.27406 |
| SCAN | 9 | 0.01189 | 0.01233 | 0.01670 | 0.03700 | 0.25504 |
| HSE06 | 9 | 0.01611 | 0.01700 | 0.02047 | 0.03500 | 0.34911 |

### a_A: literature values

| Solid | LDA | PBE | PBEsol | M06-L | SCAN | HSE06 |
| --- | --- | --- | --- | --- | --- | --- |
| AlP | 5.4340 | 5.5080 | 5.4700 | 5.4570 | 5.4750 | 5.4820 |
| BN | 3.5820 | 3.6260 | 3.6070 | 3.6030 | 3.6050 | 3.6010 |
| BP | 4.4920 | 4.5490 | 4.5210 | 4.5190 | 4.5280 | 4.5300 |
| C | 3.5320 | 3.5720 | 3.5550 | 3.5510 | 3.5510 | 3.5490 |
| LiCl | 4.9660 | 5.1510 | 5.0640 | 5.1170 | 5.1090 | 5.1070 |
| MgO | 4.1640 | 4.2540 | 4.2140 | 4.2000 | 4.1940 | 4.2080 |
| MgS | 5.1390 | 5.2310 | 5.1820 | 5.1790 | 5.1950 | 5.2070 |
| Si | 5.4050 | 5.4700 | 5.4310 | 5.4250 | 5.4330 | 5.4440 |
| SiC | 4.3300 | 4.3810 | 4.3580 | 4.3280 | 4.3530 | 4.3530 |

### B0_GPa: matched errors

Solids: AlP, BN, BP, C, LiCl, MgO, MgS, Si, SiC.

| Method | N | ME | MAE | RMSE | MaxAE | MARE (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Skala GAPW_XC-GTH | 9 | 13.10994 | 14.78472 | 18.54134 | 34.98723 | 10.39364 |
| Skala GAPW-AE | 9 | 16.55604 | 16.55604 | 20.40021 | 38.05495 | 10.41379 |
| Skala hybrid direct | 9 | 18.99265 | 20.36388 | 25.03183 | 40.12665 | 11.95837 |
| Skala hybrid +1c | 9 | 18.22250 | 19.59367 | 23.98417 | 38.05495 | 11.46466 |
| LDA | 9 | 3.12222 | 4.65556 | 6.78339 | 13.70000 | 2.95003 |
| PBE | 9 | -12.95556 | 12.95556 | 14.27048 | 23.60000 | 8.62217 |
| PBEsol | 9 | -4.72222 | 4.72222 | 6.45334 | 15.30000 | 3.52567 |
| M06-L | 9 | -1.72222 | 3.45556 | 4.14045 | 7.70000 | 2.75104 |
| SCAN | 9 | 0.43333 | 3.07778 | 3.75514 | 6.30000 | 1.80644 |
| HSE06 | 9 | 1.72222 | 5.85556 | 7.54932 | 15.50000 | 3.64562 |

### B0_GPa: literature values

| Solid | LDA | PBE | PBEsol | M06-L | SCAN | HSE06 |
| --- | --- | --- | --- | --- | --- | --- |
| AlP | 89.4000 | 82.4000 | 86.6000 | 91.1000 | 89.5000 | 94.5000 |
| BN | 402.1000 | 374.1000 | 388.0000 | 391.7000 | 394.5000 | 399.6000 |
| BP | 175.2000 | 162.0000 | 169.5000 | 170.1000 | 171.1000 | 173.0000 |
| C | 467.0000 | 434.3000 | 451.2000 | 452.8000 | 459.6000 | 468.8000 |
| LiCl | 40.8000 | 31.7000 | 35.0000 | 35.0000 | 36.3000 | 35.1000 |
| MgO | 171.7000 | 149.4000 | 157.7000 | 165.3000 | 169.6000 | 165.0000 |
| MgS | 82.8000 | 73.8000 | 78.0000 | 81.5000 | 80.1000 | 78.8000 |
| Si | 96.2000 | 88.6000 | 95.3000 | 97.7000 | 99.1000 | 97.6000 |
| SiC | 228.7000 | 212.9000 | 222.0000 | 225.1000 | 229.9000 | 228.9000 |

## Schimka2011

[Primary source](https://doi.org/10.1063/1.3524336)

Table III: lattice constants only; MgS absent. Bulk-modulus values were not extracted from its separate supplement.

### a_A: matched errors

Solids: AlN, AlP, BN, BP, C, LiCl, MgO, Si, SiC.

| Method | N | ME | MAE | RMSE | MaxAE | MARE (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Skala GAPW_XC-GTH | 9 | -0.07914 | 0.07914 | 0.08925 | 0.16391 | 1.69214 |
| Skala GAPW-AE | 9 | -0.06601 | 0.06601 | 0.07624 | 0.15650 | 1.43146 |
| Skala hybrid direct | 9 | -0.08070 | 0.08070 | 0.08878 | 0.15417 | 1.74170 |
| Skala hybrid +1c | 9 | -0.08069 | 0.08069 | 0.08877 | 0.15414 | 1.74173 |
| PBE | 9 | 0.04378 | 0.04378 | 0.04789 | 0.07600 | 0.95839 |
| PBEsol | 9 | 0.01133 | 0.01356 | 0.01640 | 0.03300 | 0.30207 |
| HSE06 | 9 | 0.01056 | 0.01322 | 0.01869 | 0.04300 | 0.27533 |
| HSEsol | 9 | -0.01122 | 0.01167 | 0.01344 | 0.02100 | 0.26695 |

### a_A: literature values

| Solid | PBE | PBEsol | HSE06 | HSEsol |
| --- | --- | --- | --- | --- |
| AlN | 4.4020 | 4.3780 | 4.3660 | 4.3510 |
| AlP | 5.5060 | 5.4720 | 5.4720 | 5.4500 |
| BN | 3.6260 | 3.6080 | 3.5980 | 3.5870 |
| BP | 4.5470 | 4.5210 | 4.5190 | 4.5040 |
| C | 3.5730 | 3.5560 | 3.5490 | 3.5380 |
| LiCl | 5.1480 | 5.0660 | 5.1150 | 5.0520 |
| MgO | 4.2600 | 4.2220 | 4.2100 | 4.1840 |
| Si | 5.4690 | 5.4360 | 5.4350 | 5.4150 |
| SiC | 4.3790 | 4.3590 | 4.3470 | 4.3340 |

## Mo2017

[Primary source](https://doi.org/10.1103/PhysRevB.95.035118)

Tables II/III: self-consistent Tao-Mo results; AlN and BN absent. Experimental references in that paper differ from Goldzak; errors here are recomputed against Goldzak.

### a_A: matched errors

Solids: AlP, BP, C, LiCl, MgO, MgS, Si, SiC.

| Method | N | ME | MAE | RMSE | MaxAE | MARE (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Skala GAPW_XC-GTH | 8 | -0.09115 | 0.09115 | 0.09881 | 0.16391 | 1.88883 |
| Skala GAPW-AE | 8 | -0.07968 | 0.07968 | 0.09043 | 0.15650 | 1.65796 |
| Skala hybrid direct | 8 | -0.09068 | 0.09068 | 0.09715 | 0.15417 | 1.88798 |
| Skala hybrid +1c | 8 | -0.09071 | 0.09071 | 0.09718 | 0.15414 | 1.88871 |
| Tao-Mo | 8 | 0.01937 | 0.01937 | 0.02158 | 0.03900 | 0.40708 |

### a_A: literature values

| Solid | Tao-Mo |
| --- | --- |
| AlP | 5.4870 |
| BP | 4.5340 |
| C | 3.5640 |
| LiCl | 5.0890 |
| MgO | 4.2090 |
| MgS | 5.1980 |
| Si | 5.4430 |
| SiC | 4.3740 |

### B0_GPa: matched errors

Solids: AlP, BP, C, LiCl, MgO, MgS, Si, SiC.

| Method | N | ME | MAE | RMSE | MaxAE | MARE (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Skala GAPW_XC-GTH | 8 | 10.38699 | 12.27112 | 15.31532 | 34.98723 | 10.57015 |
| Skala GAPW-AE | 8 | 14.57882 | 14.57882 | 18.36251 | 38.05495 | 10.67389 |
| Skala hybrid direct | 8 | 17.32001 | 18.86265 | 23.95639 | 40.12665 | 12.41154 |
| Skala hybrid +1c | 8 | 16.45359 | 17.99616 | 22.71867 | 38.05495 | 11.85611 |
| Tao-Mo | 8 | -3.31250 | 4.26250 | 5.51645 | 10.90000 | 2.53216 |

### B0_GPa: literature values

| Solid | Tao-Mo |
| --- | --- |
| AlP | 89.3000 |
| BP | 171.5000 |
| C | 442.4000 |
| LiCl | 36.2000 |
| MgO | 174.5000 |
| MgS | 79.8000 |
| Si | 97.1000 |
| SiC | 220.0000 |

## PeriodicGFN2Manuscript

Own periodic GFN manuscript, Table S3; not the original GFN2-xTB paper.

Nine-solid intersection: the manuscript LC10 excludes MgO/LiH, whereas this selected set excludes LiF/LiH. No bulk-modulus values are available in the curated table.

### a_A: matched errors

Solids: AlN, AlP, BN, BP, C, LiCl, MgS, Si, SiC.

| Method | N | ME | MAE | RMSE | MaxAE | MARE (%) |
| --- | --- | --- | --- | --- | --- | --- |
| Skala GAPW_XC-GTH | 9 | -0.08314 | 0.08314 | 0.09368 | 0.16391 | 1.72694 |
| Skala GAPW-AE | 9 | -0.07032 | 0.07032 | 0.08249 | 0.15650 | 1.46588 |
| Skala hybrid direct | 9 | -0.08307 | 0.08307 | 0.09160 | 0.15417 | 1.73801 |
| Skala hybrid +1c | 9 | -0.08305 | 0.08305 | 0.09159 | 0.15414 | 1.73764 |
| GFN1-xTB | 9 | 0.04968 | 0.13599 | 0.16141 | 0.26290 | 2.82707 |
| GFN2-xTB | 9 | -0.00987 | 0.06589 | 0.09420 | 0.24190 | 1.36417 |

### a_A: literature values

| Solid | GFN1-xTB | GFN2-xTB |
| --- | --- | --- |
| AlN | 4.6121 | 4.3128 |
| AlP | 5.5741 | 5.5450 |
| BN | 3.6252 | 3.6302 |
| BP | 4.5785 | 4.5501 |
| C | 3.5585 | 3.5549 |
| LiCl | 4.8091 | 4.8301 |
| MgS | 5.0625 | 5.1442 |
| Si | 5.6187 | 5.4462 |
| SiC | 4.5234 | 4.4127 |

## Alternate experimental B0 sensitivity

Use the three alternate values listed in Goldzak Table III simultaneously, keeping all other primary references and all ten solids unchanged. This is a reference-choice sensitivity test, not another fitted result.

| Method | N | MAE B0 primary (GPa) | MAE B0 alternate (GPa) |
| --- | --- | --- | --- |
| Skala GAPW_XC-GTH | 10 | 16.6134 | 15.6134 |
| Skala GAPW-AE | 10 | 16.8894 | 15.8894 |
| Skala hybrid direct | 10 | 21.7692 | 20.7692 |
| Skala hybrid +1c | 10 | 21.0744 | 20.0744 |
| HF | 10 | 15.2300 | 13.8900 |
| MP2 | 10 | 5.7100 | 6.8100 |
| SCS-MP2 | 10 | 4.1000 | 5.6400 |
| SOS-MP2 | 10 | 3.8800 | 5.2600 |

## Interpretation and limits

- Every native lattice constant is below the selected experimental reference. The contraction is systematic, not confined to one outlier.
- Similarity among representations does not establish accuracy against experiment or validate the implementation independently. Functional, basis/pseudopotential and numerical contributions remain entangled.
- The hybrid one-center correction changes a0 by at most 0.000542 angstrom; this does not remove the contraction. Its maximum B0 effect is 5.752 GPa (MgO).
- N=10 one-center averages include two reused AE zeros. The N=8 independent comparison is in eos-selected-method-statistics-nonshared.csv.
- Endpoint omission shifts a0 by up to 0.002444 angstrom and B0 by up to 9.997 GPa. These are fit-window sensitivity tests, not cutoff/k-point error bars.
- Goldzak lists alternate experimental B0 values for BN (410.2 versus 388.5 GPa), BP (168.0 versus 176.5) and MgO (169.8 versus 173.0). Reference sensitivity must not be confused with numerical convergence.
- Cohesive energies are not compared: no complete selected atomic-reference analysis is being claimed.

## Aggregate-only literature is not a matched ranking

Goldzak Table I lists PBE/PBEsol/SCAN MAEs of 0.061/0.030/0.030 angstrom and 12.2/7.8/7.4 GPa. These coincide with the 44-solid statistics in Tran, Stelzl and Blaha Table II; they must not be treated as newly calculated errors on this ten-solid set. The old aggregate metadata labeled them Goldzak12; this scope has been corrected to contextual/unspecified in that file.
[Tran et al., JCP 144, 204120 (2016)](https://doi.org/10.1063/1.4948636). Its system-resolved supplement was not retrieved; no missing numbers were inferred from the aggregate.
The system-resolved Zhang comparison above provides an independently sourced matched PBE/PBEsol/SCAN comparison instead.

## Reproduction

`python3 -B benchmarks/Goldzak12/scripts/compare_selected_literature.py`

Inputs: selected fits, selection JSON and five curated reference CSVs. All comparison CSVs and this report are generated offline. No CP2K jobs, energy ledgers or acceptance criteria are changed.
