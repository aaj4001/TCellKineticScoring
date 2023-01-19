# TCellKineticScoring
Allows to score gene programs in Kinetic Trajectory datasets of CD8+ T cell differentiation to contextualize T cell gene programs and infer differentiation states (see Jaiswal et al Cancer Cell 2022 for example use case)

Signature Score Calculation for Hacohen Pseudotime Clusters done on kinetic trajectories -
<ul>
  <li>LCMV Arm (Sarkar et al 2008 GSE10239/Doering et al 2012 GSE41867)</li>
  <li>LCMV cl 13 (Doering et al 2012 GSE41867)</li>
  <li>SV40TAG tumorigenesis (Philip et al 2017 GSE89307)</li>
  <li>Vaccinia virus infection (Pan et al 2017 GSE79805)</li>
</ul>

Signature Score is calculated by filtering each dataset on probes of interest, centering all gene expression across each sample (row z score) and then averaging centered expression across all genes (column mean) to compute a signature score for each sample. Signature_Up - Signature_Down <p>

Random Score is calculated by doing the same for an equal number of random probes <p>

Statistics compare Signature Score expression to Random Score expression (t test) <p>

<section>
    <h2>README v1.0</h2>
    <ul>
      <li>Moved Signature Score calculation onto one file instead of multiple for ease of reading</li>
      <li>Condensed random Score (took mean so that # of scores = # of signature Scores)</li>
      <li>Lme instead of T test to compare signature scores and random scores</li>
      <li>Added Slope analysis of Kupper Vaccinia enrichment compared to random probes</li>
    </ul>    
</section>
