# KAPA HyperPlus WGS Preparation for Illumina short-read sequencing
Generating short reads to possibly use in polishing the reference genome assembly. This prep will be added to Remy's big black sea bass WGS library which is getting submitted at the end of the month.
## Date: September 11, 2024
## Who: Annabel Hughes

**Following the same [protocol](https://remygatins.github.io/GatinsLabProtocols/lab_molec_illuminalibraryprep.html) we used in the lab to prepare the WGS libraries for the black sea bass.**

1. Dilution: need 250 ng in 11.7 uL (final concentration=21.4ng/uL)

| Sample ID | DNA Input (uL) | EB (Tris-HCL Ph=8.2) Input (uL) |
|-----------|-----------|----------------------------|
| HCI_CUR_092401_GFF | 7.86 | 3.84 |

2. Fragmentation: 37C for 06:45
3. End Repair and A-Tailing: same as protocol
4. Adapter Ligation: **Y-in line adapter #1**
5. Library Amplification: **i5 adapter #4** and **i7 adapter #1**
6. Size Selection: lower cut only (0.7x)

|8x Library (uL)| KAPA Pure Beads (uL)|
|------------|---------------------|
|16.7|11.7|

## Qubit Results:

Final concentration=22.0 ng/uL

## Tapestation Results:
Used a D1000 tape. This tapestation run looked pretty wonky... but I think it may be because I used a tape that had been open for a few days/weeks. We still sequenced it!
![plot](https://github.com/amhughes8/Hciliaris_Genome/blob/main/photos/tapestation_result_HCI_CUR_092401_WGS.png)
