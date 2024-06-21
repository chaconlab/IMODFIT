This tutorial is a practical guide for learning how to flexibly fit atomic models into low/medium resolution EM maps with IMODfit. We included some working examples with simple instructions and a benchmark setd. This tutorial is divided into three parts:</p>
<ul>
<li><a href="#Basic_fitting">Basic Flexible Fitting I. From open to close</a></li>
<li><a href="#Basic_fitting2">Basic Flexible Fitting II. From close to open</a></li>
<li><a href="#Bench52">The 52 proteins benchmark</a></li>
</ul>
<p> </p>

### Basic Flexible Fitting I. From open to close

First download, uncompress, and untar the corresponding tutorial file. You will have all the necessary files to follow the tutorial. We recommend to create a new working directory, so all the output files will be stored on it.</p>
<p>To illustrate the basic procedure and the method performance, we are going to fit the high resolution structure of GroEL monomer in an open state ( <a style="background-color: #dbe5f8;" title="Click here to view in Jmol" href="http://chaconlab.org.previewdns.com/"> <img src="images/sbg/jmol_icon.bmp" height="14" border="0" />1aon</a> , cyan) into a 10Å resolution simulated EM <a href="media/files/1oel.ccp4.gz">closed map</a> (grey) obtained from the closed structure ( <a style="background-color: #dbe5f8;" title="Click here to view in Jmol" href="http://chaconlab.org.previewdns.com/"> <img src="images/sbg/jmol_icon.bmp" height="14" border="0" />1oel</a> , yellow).</p>
<table border="0" cellspacing="2" cellpadding="0" align="center">
<tbody>
<tr>
<td><a href="images/1aonTub_354.jpg"><img style="border: 0;" title="GroEL open (1aon) (Click to enlarge)" src="images/sbg/imodfit/1aonTub_177.jpg" width="177" height="232" border="0" /> </a></td>
<td align="center" width="200">F I T T I N G<br />(open −−&gt; closed)</td>
<td><a href="images/1oelTubMapI_354.jpg"><img style="border: 0;" title="10Å resolution simulated EM map  										(from 1oel) (Click to enlarge)" src="images/sbg/imodfit/1oelTubMapI_177.jpg" width="177" border="0" /> </a></td>
</tr>
</tbody>
</table>
<p>To perform the fitting just type at the command prompt:</p>
<div class="box-note">imodfit 1aon.pdb 1oel.ccp4 10 0 -t</div>
<p>where 1aon.pdb is the initial structure, 1oel.ccp4 is the target map, 10 is the resolution in Angstroms, and 0 is the density cutoff to take into account map densities. The −t option enables the output of the PDB movie with the fitting trajectory. Here is the screen output:</p>
<pre>imodfit&gt;
imodfit&gt; Welcome to iMODFIT v1.28
imodfit&gt;
imodfit&gt; Model PDB file: 1aon.pdb
molinf&gt; Protein   1 chain  1   segment  1 residues: 524 atoms: 3847
molinf&gt; SUMMARY:
molinf&gt; Number of Molecules ... 1
molinf&gt; Number of Chain ....... 1
molinf&gt; Number of Segments .... 1
molinf&gt; Number of Groups ...... 524
molinf&gt; Number of Atoms ....... 3847
molinf&gt;
imodfit&gt; Coarse-Graining model: Full-Atom (no coarse-graining)
imodfit&gt; Selected model number of residues: 524
imodfit&gt; Selected model number of (pseudo)atoms: 3847
imodfit&gt; Target Map file: 1oel.ccp4
imodfit&gt; Best filtration method: 2 FT(x10)=0.090s Kernel(x10)=0.080s
imodfit&gt; Number of Inter-segment coords: 0 (Rot+Trans)
imodfit&gt; Number of Internal Coordinates: 1033 (Hessian rank)
imodfit&gt; Range of used modes: 1-206 (19.9%)
imodfit&gt; Number of excited/selected modes: 4(nex)
imodfit&gt;
imodfit&gt;  Iter     score     Corr. NMA   NMA_time
imodfit&gt;     0  0.336409  0.663591   0   4.33 sec
imodfit&gt;   157  0.319447  0.680553   1   4.52 sec
.................................................
imodfit&gt;  4219  0.024751  0.975249  17   4.15 sec
imodfit&gt; 10000  0.017593  0.982407            END
imodfit&gt;
imodfit&gt; Movie file:                              imodfit_movie.pdb
imodfit&gt; Final Model:                            imodfit_fitted.pdb
imodfit&gt; Score file:                              imodfit_score.txt
imodfit&gt; Log file:                                      imodfit.log
imodfit&gt;
imodfit&gt; Success! Time elapsed  00h. 04' 04''
imodfit&gt; Bye!
</pre>
<p>The flexibly fitted structure is: <a style="background-color: #dbe5f8;" title="Click here to view in Jmol" href="http://chaconlab.org.previewdns.com/"> <img src="images/sbg/jmol_icon.bmp" height="14" border="0" /> <b>imodfit_fitted.pdb</b></a>.</p>
<p>iMODFIT also outputs the following files: </p>
<ul>
<li><i>imodfit_movie.pdb </i>--&gt; fitting trajectory</li>
<li><i>imodfit_score.txt </i>--&gt; score file to check for convergence</li>
<li><i>imodfit.log</i>--&gt; used command log</li>
</ul>
<p>Below some fitting trajectory snapshots (cyan) are represented simultaneously with the target structure (yellow). The final snapshot with the fitted structure is shown on the right.</p>
<table style="width: 400px;" border="0" cellspacing="2" cellpadding="0" align="center">
<tbody>
<tr>
<td><a href="images/sbg/imodfit/1aonTubMapFTub_initial_354.jpg"><img style="border: 0;" title="Initial Model (Click to enlarge)" src="images/sbg/imodfit/1aonTubMapFTub_initial_177.jpg" height="116" border="0" /></a></td>
<td><a href="images/sbg/imodfit/1aonTubMapFTub_05_354.jpg"><img style="border: 0;" title="Intermediate Model (Click to enlarge)" src="images/sbg/imodfit/1aonTubMapFTub_05_177.jpg" height="116" border="0" /></a></td>
<td><a href="images/sbg/imodfit/1aonTubMapFTub_18_354.jpg"><img style="border: 0;" title="Intermediate Model (Click to enlarge)" src="images/sbg/imodfit/1aonTubMapFTub_18_177.jpg" height="116" border="0" /></a></td>
<td><a href="images/sbg/imodfit/1aonTubMapFTub_31_354.jpg"><img style="border: 0;" title="Intermediate Model (Click to enlarge)" src="images/sbg/imodfit/1aonTubMapFTub_31_177.jpg" height="116" border="0" /></a></td>
<td><a href="images/sbg/imodfit/1aonTubMapFTub_40_354.jpg"><img style="border: 0;" title="Intermediate Model (Click to enlarge)" src="images/sbg/imodfit/1aonTubMapFTub_40_177.jpg" height="116" border="0" /></a></td>
<td><a href="images/sbg/imodfit/1aonTubMapFTub_final_354.jpg"><img style="border: 0;" title="Final Fitted Model (Click to enlarge)" src="images/sbg/imodfit/1aonTubMapFTub_final_177.jpg" height="116" border="0" /></a></td>
</tr>
</tbody>
</table>
<p>The fitting result is only 1.78Å Cα RMSD from the target structure, and the final correlation was high: 0.982. The quality of fitness and the excellent secondary structure maintenance can be appreciated in the flash movies below (front and rear views in left and right, respectively). Note that in none case any secondary structure constraint was used.</p>
<table border="0" cellspacing="0" cellpadding="0" align="center">
<tbody>
<tr>
<td><video controls="controls" width="320" height="259">
  <source src="images/video/imodfit/1aonCart_1oelMapFTub.mp4" type="video/mp4" />
  Your browser does not support HTML5 video.
  </video></td>
<td><video controls="controls" width="320" height="259">
  <source src="images/video/imodfit/1aonCart_1oelMapFTub2.mp4" type="video/mp4" />
  Your browser does not support HTML5 video.
  </video></td>
</tr>
</tbody>
</table>
<p>For visualizing the results use your favorite program. To play the trajectory movie ( <i>imodfit_movie.pdb</i> ) we recomend <a href="http://www.ks.uiuc.edu/Research/vmd/">VMD</a>; but you can see it in <a style="background-color: #dbe5f8;" title="Click here to view in Jmol" href="http://chaconlab.org.previewdns.com/"> <img src="images/sbg/jmol_icon.bmp" height="14" border="0" /> Jmol</a>. The images and the movie were created using VMD with the POVray rendering engine.</p>
<p>Once you have run iMODFIT you should check for convergence. To this end just plot the score (i.e. 1−correlation) as a function of the iteration number using GNUplot with the <i>imodfit_score.txt</i> file:</p>
<pre>&gt; gnuplot -persist
&gt; gnuplot&gt; plot "imodfit_score.txt" u 1:2 w l
&gt; gnuplot&gt; exit
</pre>
<table border="0" cellspacing="2" cellpadding="0" align="center">
<tbody>
<tr>
<td><a href="images/imodfit_score_640.jpg"><img style="border: 0;" title="Score profile  											(Click to Enlarge)" src="images/sbg/imodfit/imodfit_score_320.jpg" width="320" border="0" /></a></td>
</tr>
</tbody>
</table>
<p>Note(1): If the slope is not approximately horizontal at 5000−10000 you should run iMODFIT again to increase the number of maximum iterations (−i option). Alternatively, you can continue the fitting process introducing the final fitted structure as the initial model.</p>
<p>Note(2): The conformational refinement is a stochastic process; thus you should not expect to obtain the same screen output and the same results between different runs.</p>
<p> </p>

### Basic Flexible Fitting II. From close to open

iMODFIT allows also the fitting from the closed structures as well. To perform this fitting just try the following command:</p>
```
imodfit 1oel.pdb 1aon.ccp4 10 0 -t -o imodfit2
```
<p>Front and rear views of the fitting results are shown in the flash movies below. The initial closed structure, 1oel, is shown in cyan, the target open one, 1aon, in yellow, and its corresponding <a href="media/files/1aon.ccp4.gz">open map</a> in grey. The −o option is added to avoid overwriting previous results.</p>
<table style="width: 400px;" border="0" cellspacing="2" cellpadding="0" align="center">
<tbody>
<tr>
<td><video controls="controls" width="320" height="259">
  <source src="images/video/imodfit/1oelCart_1aonMapFTub.mp4" type="video/mp4" />
  Your browser does not support HTML5 video.
  </video></td>
<td><video controls="controls" width="320" height="259">
  <source src="images/video/imodfit/1oelCart_1aonMapFTub2.mp4" type="video/mp4" />
  Your browser does not support HTML5 video.
  </video></td>
</tr>
</tbody>
</table>
<p>The final structure is again very close to the target, i.e. 1.74Å Cα RMSD (corr.=0.985). As in above case, the quality of fitness and the excellent Secondary Structure maintenance are evident. Also, you can observe the trajectory interactively: <a style="background-color: #dbe5f8;" title="Click here to view in Jmol" href="http://chaconlab.org.previewdns.com/"> <img src="images/sbg/jmol_icon.bmp" height="14" border="0" /> Jmol</a>.</p>
<p> </p>

### The 52 proteins benchmark

<p>To test our methodology we build a benchmark formed by 52 simulated fitting problems, comprising a wide variety of macromolecular motions. To this end, we downloaded from the molecular motions database <a href="http://molmovdb.org/cgi-bin/browse.cgi">MolMovDB</a> 26 different protein pairs with displacements greater than 2Å Cα RMSD and sizes ranging from 100 to 1000 amino acids. Only those structures with less than 3% Ramachandran outliers (<a href="http://molprobity.biochem.duke.edu">Molprobity</a>) and without broken chains and missing atoms were admitted. The average displacement was 6.9Å with a standard deviation of 3.9Å. Note that each open/closed protein pair represents two different fitting problems, i.e. from open to closed and vice versa. It can be downloaded either from the iMODFIT "Install" page or directly from <a class="dwnld" href="downloads/Tools/iMODFIT/bench52.tar.gz">here</a>.</p>
<table class="text" border="0" cellspacing="0" cellpadding="6" align="center">
<tbody>
<tr>
<td>
<table class="text" border="0" cellspacing="0" cellpadding="3" align="center" bgcolor="#f5f5f5">
<tbody>
<tr>
<td><i>Open</i></td>
<td><i>Closed</i></td>
<td><i>Motion</i></td>
<td><i>Name</i></td>
<td> </td>
<td><i>#Residues</i></td>
<td><i>Cα RMSD (Å)</i></td>
</tr>
<tr>
<td><b>1l5e</b></td>
<td><b>1l5b</b></td>
<td><i>[D-h-2]</i></td>
<td>Cyanovirin-N</td>
<td> </td>
<td align="center">101</td>
<td align="center">8.85</td>
</tr>
<tr>
<td><b>1e7xA</b></td>
<td><b>1dzsB</b></td>
<td><i>[F-?-2]</i></td>
<td>Virus MS2 coat protein</td>
<td> </td>
<td align="center">129</td>
<td align="center">3.57</td>
</tr>
<tr>
<td><b>1cfd</b></td>
<td><b>1cfc</b></td>
<td><i>[D-h-2]</i></td>
<td>Calmodulin</td>
<td> </td>
<td align="center">148</td>
<td align="center">10.22</td>
</tr>
<tr>
<td><b>1cbuB</b></td>
<td><b>1c9kB</b></td>
<td><i>[F-s-2]</i></td>
<td>Adenosylcobinamide Kinase</td>
<td> </td>
<td align="center">180</td>
<td align="center">3.52</td>
</tr>
<tr>
<td><b>1ex6</b></td>
<td><b>1ex7</b></td>
<td><i>[D-h-2]</i></td>
<td>Guanylate Kinase (GDK)</td>
<td> </td>
<td align="center">186</td>
<td align="center">4.58</td>
</tr>
<tr>
<td><b>4ake</b></td>
<td><b>1ake</b></td>
<td><i>[D-h-2]</i></td>
<td>Adenylate Kinase (ADK)</td>
<td> </td>
<td align="center">214</td>
<td align="center">8.26</td>
</tr>
<tr>
<td><b>1gggA</b></td>
<td><b>1wdnA</b></td>
<td><i>[D-h-2]</i></td>
<td>Glutamine Binding Protein</td>
<td> </td>
<td align="center">220</td>
<td align="center">5.34</td>
</tr>
<tr>
<td><b>2lao</b></td>
<td><b>1lst</b></td>
<td><i>[D-h-2]</i></td>
<td>Lysine/Arginine/Ornithine (LAO) binding protein</td>
<td> </td>
<td align="center">238</td>
<td align="center">8.67</td>
</tr>
<tr>
<td><b>1urp</b></td>
<td><b>2dri</b></td>
<td><i>[D-h-2]</i></td>
<td>Ribose Binding Protein</td>
<td> </td>
<td align="center">271</td>
<td align="center">7.69</td>
</tr>
<tr>
<td><b>1ram</b></td>
<td><b>1leiA</b></td>
<td><i>[D-?-2]</i></td>
<td>NF-kappa B</td>
<td> </td>
<td align="center">273</td>
<td align="center">4.96</td>
</tr>
<tr>
<td><b>5at1</b></td>
<td><b>8atc</b></td>
<td><i>[S-a-2]</i></td>
<td>Aspartate Carbamoyltransferase</td>
<td> </td>
<td align="center">310</td>
<td align="center">2.41</td>
</tr>
<tr>
<td><b>1ckmA</b></td>
<td><b>1ckmB</b></td>
<td><i>[D-h-2]</i></td>
<td>mRNA capping enzyme</td>
<td> </td>
<td align="center">317</td>
<td align="center">4.35</td>
</tr>
<tr>
<td><b>3dap</b></td>
<td><b>1dap</b></td>
<td><i>[D-h-2]</i></td>
<td>Diaminopimelic Acid Dehydrogenase</td>
<td> </td>
<td align="center">320</td>
<td align="center">5.81</td>
</tr>
<tr>
<td><b>1bp5</b></td>
<td><b>1a8e</b></td>
<td><i>[D-h-2]</i></td>
<td>Transferrins (N-terminal lobe)</td>
<td> </td>
<td align="center">329</td>
<td align="center">12.16</td>
</tr>
<tr>
<td><b>1jqj</b></td>
<td><b>2pol</b></td>
<td><i>[D-s-2]</i></td>
<td>Beta DNA polymerase III</td>
<td> </td>
<td align="center">366</td>
<td align="center">2.81</td>
</tr>
<tr>
<td><b>1omp</b></td>
<td><b>1anf</b></td>
<td><i>[D-h-2]</i></td>
<td>Maltodextrin Binding Protein (MBP)</td>
<td> </td>
<td align="center">370</td>
<td align="center">7.23</td>
</tr>
<tr>
<td><b>8adh</b></td>
<td><b>6adh</b></td>
<td><i>[D-s-2]</i></td>
<td>Alcohol Dehydrogenase (ADH)</td>
<td> </td>
<td align="center">374</td>
<td align="center">2.43</td>
</tr>
<tr>
<td><b>9aat</b></td>
<td><b>1ama</b></td>
<td><i>[D-s-2]</i></td>
<td>Aspartate Amino Transferase (AAT)</td>
<td> </td>
<td align="center">401</td>
<td align="center">2.21</td>
</tr>
<tr>
<td><b>1bnc</b></td>
<td><b>1dv2</b></td>
<td><i>[D-h-2]</i></td>
<td>Biotin carboxylase</td>
<td> </td>
<td align="center">452</td>
<td align="center">5.38</td>
</tr>
<tr>
<td><b>1rkm</b></td>
<td><b>2rkm</b></td>
<td><i>[D-h-2]</i></td>
<td>Oligopeptide-binding protein</td>
<td> </td>
<td align="center">517</td>
<td align="center">5.77</td>
</tr>
<tr>
<td><b>1sx4</b></td>
<td><b>1oel</b></td>
<td><i>[C----]</i></td>
<td>GroEL</td>
<td> </td>
<td align="center">524</td>
<td align="center">15.83</td>
</tr>
<tr>
<td><b>1i7d</b></td>
<td><b>1d6m</b></td>
<td><i>[D-?-2]</i></td>
<td>DNA Topoisomerase III</td>
<td> </td>
<td align="center">620</td>
<td align="center">4.16</td>
</tr>
<tr>
<td><b>8ohm</b></td>
<td><b>1cu1</b></td>
<td><i>[D-h-2]</i></td>
<td>HCV Helicase</td>
<td> </td>
<td align="center">645</td>
<td align="center">6.06</td>
</tr>
<tr>
<td><b>1lfg</b></td>
<td><b>1lfh</b></td>
<td><i>[D-h-2]</i></td>
<td>Lactoferrin (LF)</td>
<td> </td>
<td align="center">691</td>
<td align="center">8.08</td>
</tr>
<tr>
<td><b>1ih7</b></td>
<td><b>1ig9</b></td>
<td><i>[C----]</i></td>
<td>Phage Rb69 DNA Polymerase</td>
<td> </td>
<td align="center">903</td>
<td align="center">7.16</td>
</tr>
<tr>
<td><b>1su4</b></td>
<td><b>1t5s</b></td>
<td><i>[C----]</i></td>
<td>Calcium ATPase pump</td>
<td> </td>
<td align="center">994</td>
<td align="center">17.99</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>
<p>MolMovDB's Motion Types:</p>
<ol style="list-style-type: lower-roman;">
<li>Motions of Fragments Smaller than Domains
<ul>
<li><em>[F-s-2] </em>--&gt; I.A. Motion is predominantly shear</li>
<li><i>[F-h-2] </i>--&gt; I.B. Motion is predominantly hinge</li>
<li><i>[F-?-2] </i>--&gt; I.C. Motion can not be fully classified at present</li>
<li><i>[F-n-2] </i>--&gt; I.D. Motion is not hinge or shear</li>
</ul>
</li>
<li>Domain Motions
<ul>
<li><i>[D-s-2] </i>--&gt; II.A. Motion is predominantly shear</li>
<li><i>[D-h-2] </i>--&gt; II.B. Motion is predominantly hinge</li>
<li><i>[D-?-2] </i>--&gt; II.C. Motion can not be fully classified at present</li>
<li><i>[D-n-2] </i>--&gt; II.D. Motion is not hinge or shear</li>
<li><i>[D-f-2] </i>--&gt; II.E. Motion involves partial refolding of tertiary structure</li>
</ul>
</li>
<li>Larger Movements than Domain Motions involving the Movement of Subunits
<ul>
<li><i>[S-a-2] </i>--&gt; III.A. Motion involves an allosteric transition</li>
<li><i>[S-n-2] </i>--&gt; III.B. Motion does not involves an allosteric transition</li>
<li><i>[C----] </i>--&gt; C----. Complex Protein Motions</li>
</ul>
</li>
</ol>
