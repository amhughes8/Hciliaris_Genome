# Oxford Nanopore Technologies Ligation Sequencing Library Preparation + Sequencing
## [Ligation sequencing DNA V14 (SQK-LSK114)](https://nanoporetech.com/document/genomic-dna-by-ligation-sqk-lsk114)
Annabel M. Hughes

March 24, 2025

### Step 1: Prepare your input DNA
I wanted to use 3 μg of DNA input in order to prepare a library sufficient for 3 loads on the PromethION. Based on my Nanodrop DNA concentration, I calculated the following:

| Sample ID | Concentration (ng/μL) | Volume to get to 3 μg |
| --------|-----------|------------|
| GFF | 46.9 | **63.97** |

### Step 2: DNA repair and end-prep
**Prep for this step:** 
1. Thaw all reagents on ice (see list below) and spin them down before first use of the day (_Note: DO NOT vortex the FFPE DNA Repair Mix or the Enzyme_)
2. Prepare 750 μL of fresh 80% ethanol in nuclease-free water

- [ ] In a 0.2 mL PCR tube, mix the following:
| Reagent | Volume |
|--------|--------|
| Input DNA | 63.97 μL |
| NEBNext FFPE DNA Repair Buffer | 3.5 μL |
| NEBNext FFPE DNA Repair Mix | 2 μL |
| Ultra II End-prep Enzyme | 3 μL |
| Ultra II End-prep reaction buffer | 3.5 μL |
|**Total**|**75.97 μL**|

- [ ] Incubate at 20°C for 30 minutes and 65°C for 30 minutes
- [ ] Resuspend AMPure XP Beads (AXP) by vortexing
- [ ] Transfer the DNA sample to a clean 1.5 mL Eppendorf DNA LoBind tube
- [ ] Add 76 μL (same volume as above) of resuspended AMPure XP Beads to the end-prep reaction and mix by flicking the tube
- [ ] Rotate by hand for 10 minutes at room temperature
- [ ] Spin down the sample and pellet on a magnet until supernatant is clear and colorless. Keep the tube on the magnet and pipette off supernatant.
- [ ] Keep the tube on the magnet and wash the beads with 300 μL of freshly prepared 80% ethanol without disturbing the pellet. Remove ethanol using a pipette and discard.
- [ ] Repeat the previous step
- [ ] Spin down and place the tube back on the magnet. Pipette off any residual ethanol. Allow to dry for ~30 seconds, but do not dry the pellet to the point of cracking.
- [ ] Remove the tube from the magnetic rack and resuspend the pellet in 62 μL (60 + 2 μL for Qubit) nuclease-free water. Incubate for 30 minutes at 37°C.
- [ ] Pellet the beads on a magnet until the eluate is clear and colorless, for at least 1 minute
- [ ] Remove and retain 62 μL of eluate into a clean 1.5 mL Eppendorf DNA LoBind tube.

**Checkpoint: Quantify 2 μL of eluted sample using a Qubit fluorometer.**

Qubit concentration: **37.4 ng/μL**

### Step 3: Adapter ligation and clean-up
**Prep for this step:** 
1. Spin down the Ligation Adapter (LA) and NEBNext Quick T4 DNA Ligase and place on ice
2. Thaw Ligation Buffer (LNB) at room temperature, spin down and mix by pipetting. Due to viscosity, vortexing this buffer is ineffective. Place on ice immediately after thawing and mixing.
3. Thaw Elution Buffer (EB) at room temperature and mix by vortexing. Then spin down and place on ice
4. Thaw Long Fragment Buffer (LFB) at room temperature and mix by vortexing. Then spin down and keep at room temperature.

- [ ] Mix the following in a 1.5 mL Eppendorf DNA LoBind tube in the following order:
| Reagent | Volume |
|--------|--------|
| DNA sample from previous step | 60 μL |
| Ligation Adapter (LA) | 5 μL |
| Ligation Buffer (LNB) | 25 μL |
| NEBNext Quick T4 DNA ligase | 10 μL |
|**Total**|**100 μL**|

- [ ] Thoroughly mix the reaction by gently flicking the tube and spin down
- [ ] Incubate the reaction for 60 minutes at room temperature
- [ ] Resuspend the AMPure XP Beads (AXP) by vortexing
- [ ] Add 40 μL of resuspended AMPure XP Beads to the reaction and mix by flicking the tube
- [ ] Rotate by hand for 5 minutes
- [ ] Spin down the sample and pellet on a magnet. Keep the tube on the magnet and pipette off the supernatant when clear and colorless
- [ ] Wash the beads by adding 250 μL Long Fragment Buffer. Spin down, then return the tube to the magnetic rack and allow the beads to pellet. Remove the supernatant using a pipette and discard.
- [ ] Repeat the previous step
- [ ] Spin down and place the tube back on the magnet. Pipette off any residual supernatant. Allow to dry for ~30 seconds but do not dry the pellet to the point of cracking
- [ ] Remove the tube from the magnetic rak and resuspend the pellet in 27 μL (25 + 2 μL for Qubit) of Elution Buffer. Spin down and incubate for 30 minutes at 37°C.
- [ ] Pellet the beads on a magnet until the eluate is clear and colorless, for at least one minute.
- [ ] Remove and retain 27 μL of eluate containing the DNA library into a clean 1.5 mL Eppendorf DNA LoBind tube. 

**Checkpoint: Quantify 2 μL of eluted sample using a Qubit fluorometer.**

Qubit concentration: **47.5 ng/μL**

### Step 4: Divide library into separate loads
Since we are planning to load and reload the PromethION flow cell 3 times, we will divide our final library into 3 separate aliquots to run individually. Washing the reloading the flow cell allows us to unclog some pores and increase the number of pores that are able to continue sequencing, thereby increasing our overall yield for each flow cell. The two subsequent washes and reloads will take place 24 and 48 hours after the initial load, respectively.

Final library volume: 25 μL
Final library concentration: 47.5 ng/μL

25 ÷ 3 = 8.3 μL
8.3 x 47.5 = 395 ng total library per load -- this is GREAT! Ultimately, I decided to split the loads up relatively evenly, with the first load acquiring the greatest volume since it will have the greatest number of pores active for sequencing.

Load 1: 9 μL --> **427.5 ng**
Load 2: 8 μL --> **380 ng**
Load 3: 8 μL --> **380 ng**

### Step 5: Priming and loading the PromethION flow cell
**Prep for this step:** 
1. Thaw the Sequencing Buffer (SB), Library Beads (LIB), Flow Cell Tether (FCT), and Flow Cell Flush (FCF) at room temperature before mixing by vortexing. Then spin down and store on ice.

**We will be preparing 3 tubes to bring with us for sequencing:**
1 - Flow cell priming mix (Flow Cell Flush + Flow Cell Tether)
2 - Library loading prep mix (Sequencing Buffer + Library Beads)
3 - DNA library 

## To prepare Tube 1: Flow cell priming mix (Flow Cell Flush + Flow Cell Tether)
| Reagent | Volume |
|--------|--------|
| Flow Cell Flush | 1170 μL |
| Flow Cell Tether | 30 μL |

## To prepare Tube 2: Library loading prep mix
| Reagent | Volume |
|--------|--------|
| Sequencing Buffer | 100 μL |
| Library Beads | 68 μL |

## To prepare Tube 3: DNA library 
I took 9 μL of my final library and left the remaining 16 μL in the 4°C fridge. This will eventually be added to the library loading prep mix, but keep it separate until just before loading!

Physical priming and loading of the flow cell was done according to manufacturer's instructions. In the end, our first flow cell ended up running for the full 72 hours it was allotted before getting restarted for a few hours with a small amount of remaining pores. We then loaded the remaining 16 μL of library onto a second flow cell after size selecting for larger fragments.






