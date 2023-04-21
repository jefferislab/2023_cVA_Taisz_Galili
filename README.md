# 2023_cVA_Taisz_Galili
Public data and code for connectomic analyses in Taisz, Galili et al. (2023)

To analyse connectomic data we relied on the `natverse` toolbox from our group: [Bates et al. (2020)](https://doi.org/10.7554/eLife.53350).
https://natverse.org/

Datasets:
1. FAFB – [Zheng et al. (2018)](https://doi.org/10.1016/j.cell.2018.06.019), [Li et al. (2020)](https://doi.org/10.1101/605634), [Dorkenwald et al. (2022)](https://doi.org/10.1038/s41592-021-01330-0)
Manual tracings of FAFB neurons were deposited at https://catmaid-fafb.virtualflybrain.org/

We used the `rcatmaid` package to query data from CATMAID, and the `fafbseg` package to query data from FlyWire.

https://natverse.org/rcatmaid/

https://github.com/natverse/fafbseg

https://flywire.ai/

2. hemibrain – [Scheffer et al. (2020)](https://doi.org/10.7554/eLife.57443)
The data is publically available at https://neuprint.janelia.org (hemibrain:v1.2.1).

We used the `neuprintr` package to query hemibrain data.

https://natverse.org/neuprintr/

Neuron identifiers (Skeleton IDs or bodyIDs) and the number of synaptic connections across cell types from both datasets can be found in Supplemental Tables 3 and 4. Our group identified neurons in the hemibrain prior to publication and contributed the annotation of all ORNs, PNs, and LH neurons in neuprint [(Schlegel, Bates et al. 2021)](https://doi.org/10.7554/eLife.66018).
Supplemental Table 5 shows specific analyses about third-order cell types of the cVA circuit related to Figure S6A.

All other data and analyses are available upon request. Lead contact: Greg Jefferis, https://github.com/jefferis

For generating image stacks from EM neurons for MIP search see the skeleton-to-MIP repository.
https://github.com/jefferislab/skeleton-to-MIP
