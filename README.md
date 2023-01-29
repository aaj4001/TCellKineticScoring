# TCellKineticScoring
Allows to score gene programs in Kinetic Trajectory datasets of CD8+ T cell differentiation to contextualize T cell gene programs and infer differentiation states (see Jaiswal et al Cancer Cell 2022 for example use case)

Signature Score Calculation for Hacohen Pseudotime Clusters done on kinetic trajectories -
<ul>
  <li>LCMV Arm (Sarkar et al 2008 GSE10239/Doering et al 2012 GSE41867)</li>
  <li>LCMV cl 13 (Doering et al 2012 GSE41867)</li>
  <li>SV40TAG tumorigenesis (Philip et al 2017 GSE89307)</li>
  <li>Vaccinia virus infection (Pan et al 2017 GSE79805)</li>
  <li><b>Human Yellow Fever vaccination (Akondy et al 2017) T cell activation microarray (GSE26347), T cell activation and memory RNAseq (GSE100745)</b></li>

</ul>

Signature Score is calculated by filtering each dataset on probes of interest, centering all gene expression across each sample (row z score) and then averaging centered expression across all genes (column mean) to compute a signature score for each sample. Signature_Up - Signature_Down <p>

Random Score is calculated by doing the same for an equal number of random probes <p>

  Statistics compare Signature Score expression to expression at preceding timepoint (lme for all timecourses except YF RNAseq, where a t-test was used). For vaccinia virus skin scarification, statistics of expression change over time calculated using lm. <p>

<section>
    <h2>README v5.0 <b>(This version was the final version used in Jaiswal et al Cancer Cell 2022)</b> </h2>
    <ul>
      <li>Added Human Yellow Fever datasets to kinetic set analysis </li>
    </ul>
    Scoring Functions: 
    <ul>
      <li>SignatureScorePlot_TRMSlopeEval - Signature Scoring, computes slope from d0-90 and p value of enrichment in Kupper TRM set using lm </li>
    </ul>
</section>
