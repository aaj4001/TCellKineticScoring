# TCellKineticScoring
Allows to score gene programs in Kinetic Trajectory datasets of CD8+ T cell differentiation to contextualize T cell gene programs and infer differentiation states (see Jaiswal et al Cancer Cell 2022 for example use case)

Signature Score Calculation for Hacohen Pseudotime Clusters done on kinetic trajectories -
<ul>
  <li>LCMV Arm (Sarkar et al 2008 GSE10239/Doering et al 2012 GSE41867)</li>
  <li>LCMV cl 13 (Doering et al 2012 GSE41867)</li>
  <li>SV40TAG tumorigenesis (Philip et al 2017 GSE89307)</li>
  <li>Vaccinia virus infection (Pan et al 2017 GSE79805)</li>
  <li>Human Yellow Fever vaccination (Akondy et al 2017) T cell activation microarray (GSE26347), T cell activation and memory RNAseq (GSE100745)</li>

</ul>

Signature Score is calculated by filtering each dataset on probes of interest, centering all gene expression across each sample (row z score) and then averaging centered expression across all genes (column mean) to compute a signature score for each sample. Signature_Up - Signature_Down <p>

Random Score is calculated by doing the same for an equal number of random probes <p>

  Statistics compare Signature Score expression to expression at preceding timepoint (lme for all timecourses except YF RNAseq, where a t-test was used). For vaccinia virus skin scarification, statistics of expression change over time calculated using lm. <p>

<section>
    <h2>README v6.0 </h2>
    <ul>
      <li>Redid kinetic derivation for microarray datasets. For probes specific for the same gene, the probe with the highest median absolute deviation (MAD) was used.</li>
      <li>Additionally, any microarray probes that have an sd of 0 (no variance) are also removed. </li>
      <li>Updated Kinetic scoring function. Quantification/Statistic function moved outside plotting function for modularity.</li>
      <li>Additionally, instead of printing p values to console, function now allows user to decide whether to print P values on output graph or not. </li>
      <li>Finally, user can choose whether to keep signature score plots on the same scale or not. </li>
    </ul>
    <b>Scoring Function: SignatureScorePlot</b>
    <ul>
      <li>KineticSets: Kinetic datasets scored (default should be KineticLists,see KineticSetDerivation file for details on how it was made)</li>
      <li>SigList: Signatures used for scoring (can be a list with sublist "Up" and "Down" signatures, or a list with just "Up" signatures)</li>
      <li>pal: Color palette used for each signature</li>
      <li>nperm: number of permutations for random score calculation</li>
      <li>SameScale: logical whether user wants all score plots on the same scale</li>
      <li>PValsAnnot: logical whether user wants to plot p values or not</li>
      <li>PValsSize: numeric - font of P values</li>
    </ul>
</section>
