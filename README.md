# NK cell cytotoxicity shapes the clonal evolution of B cell leukaemia 

**Link:** XXX

**Authors:**\
Michelle C. Buri, Mohamed R. Shoeb, Aleksandr Bykov, Peter Repiscak, Hayeon Baik, Alma Dupanovic, Faith O. David, Boris Kovacic, Faith Hall-Glenn, Sara Dopa, Jos Urbanus, Lisa Sippl, Susanne Stofner, Dominik Emminger, Jason Cosgrove, Dagmar Schinnerl, Anna R. Poetsch, Manfred Lehner, Xaver Koenig, Leïla Perié, Ton N. Schumacher, Dagmar Gotthardt, Florian Halbritter, Eva M. Putz

## Abstract
The term cancer immunoediting describes the dual role by which the immune system can suppress and promote tumour growth and is divided into three phases: elimination, equilibrium and escape. The role of NK cells has mainly been attributed to the elimination phase. Here we show that NK cells play a role in all three phases of cancer immunoediting. Extended co-culturing of DNA barcoded mouse BCR/ABLp185+ B acute lymphoblastic leukaemia cells with NK cells allowed for a quantitative measure of NK cell-mediated immunoediting. Whereas most tumour cell clones were efficiently eliminated by NK cells, a certain fraction of tumour cells harboured an intrinsic primary resistance. Furthermore, DNA barcoding revealed tumour cell clones with secondary resistance, which stochastically acquired resistance to NK cells. NK cell cytotoxicity put a selective pressure on B-ALL cells inducing primary and secondary resistance, while resistant tumour cells were characterised by a full-blown IFN- signature. Besides well-known regulators of immune evasion, our analysis of NK resistant tumour cells revealed the upregulation of novel genes, including Ly6a, which we found to drive NK cell resistance in leukaemic cells. We further translated our findings to the human system and showed that high LY6E expression on tumour cells impaired the physical interaction with NK cells and led to worse prognosis in leukaemia. Our results demonstrate that tumour cells are actively edited by NK cells during the equilibrium phase and use different avenues to escape NK cell-mediated eradication.
________

### Usage of the code:
This is a supplementarry code that was used to analyze the data produced in the current research project. This code relies on usage of dockerized RStudio server. Build the image using the `_set_up_the_environment.sh` and run the `dockRstudio_v4.2.0_dev_run.sh`. You can connect to the created RStudioServer by typing the IP address of the HPC you are using followed by the port number (48900 in this case) in your browser, as follows: `localhost:48900`.

The ATAC-seq is placed in a separate [repository](https://github.com/TumorImmunoEditingLab/NK_tumor_resistance_ATAC-seq).
The code was developed and tested on Linux Debian v8 (jessie).


