# Qiagen Blood and Tissue DNA extractions for Illumina + ONT
## Date: September 5, 2024
## Who: Annabel Hughes

1. Using four tissue types to see which gives us the best quality extraction:

| Tissue Type | Preservation Method | Freezer Storage | Sample ID |
|-----------|-------------------|----------------|----------|
| fin | ethanol | -80C | HCI_CUR_092401_FIE |
| fin | flash frozen | -80C | HCI_CUR_092401_FFF |
| gill | flash frozen | -80C | HCI_CUR_092401_GFF |
| skin | ethanol | -80C | HCI_CUR_092401_SIE |

2. Incubated for 1 hour at 56C (vortexed and spun every 20 mins).
3. Eluted in 200 uL EB (Tris-HCL pH=8.2). The EB was not pre-heated.
4. Repeated elution with original 200 uL

## Qubit Results:
| Sample ID | Concentration [ng/uL] |
|---------|---------------------|
| HCI_CUR_092401_FIE | 17.9 |
| HCI_CUR_092401_FFF | 9.05 |
| HCI_CUR_092401_GFF | 31.8 |
| HCI_CUR_092401_SIE | 4.88 |

## Nanodrop Results:
| Sample ID | Concentration [ng/uL] | 260/280 | 260/230 |
|---------|--------------------|---------|---------|
| HCI_CUR_092401_FIE | 32.5 | 2.20 | 1.79 |
| HCI_CUR_092401_FFF | 19.3 | 1.87 | 1.38 |
| HCI_CUR_092401_GFF | 46.9 | 1.95 | 1.90 |
| HCI_CUR_092401_SIE | 12 | 2.15 | 1.05 |

## Gel Results:
2% gel, 90 V for 35 mins
![plot](photos/gel_results.png)

## Tapestation Results:
Ran a gDNA tape with FIE, FFF, and GFF. GFF looked the past through all quality checks:
![plot](photos/tapestation_results_gDNA.png)

