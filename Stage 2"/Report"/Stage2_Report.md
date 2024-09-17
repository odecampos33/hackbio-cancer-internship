**Gene Expression Patterns in Glioblastoma: Heatmap Visualization and Pathway Enrichment Analysis**

@Hackbio @slack 

 _Authors:_

- Chaimaa BenMohamed
- Charlotte Chinwendu Iwuji
- Chigozie Nkwocha
- Igwebuike Vivian
- Opeyemi De Campos
- Reem  Atawia****

**Overview**

This task focuses on visualizing and interpreting a gene expression dataset related to glioblastoma, the most common primary malignant brain tumor. The aim is to generate heatmaps, identify differentially expressed genes (DEGs), and perform functional enrichment analysis to interpret biological processes that contribute to the disease.

**Protein with unfavorable prognosis in glioblastoma,**

Numerous genes are implicated in the pathogenesis of glioblastoma, including the EGF (Epidermal Growth Factor) gene, which encodes the Epidermal Growth Factor Receptor (EGFR) protein that plays a critical role in cell signaling pathways regulating cell division and survival (Brosseau et al., 2015; Huang et al., 2023; Chen et al., 2016; Singh et al., 2023).

![Figure 1](https://github.com/odecampos33/hackbio-cancer-internship/blob/main/Stage%202%22/Images%22/EGFR_signaling.png?raw=true)
Figure 1: EGFR signaling pathway 
**Data Preprocessing and Visualization**

Using the glioblastoma dataset, containing over 500 differentially expressed genes from 10 patients, the first step involved preprocessing, which included data normalization and exploration. Heatmaps were generated using a diverging and sequential color palette to improve interpretation. The diverging palette effectively visualized the difference between upregulated and downregulated genes, while the sequential palette allowed a smoother transition in gene expression intensities. Both palettes were essential for making high, low, and medium values easily distinguishable.

![Figure 2](https://github.com/odecampos33/hackbio-cancer-internship/blob/main/Stage%202%22/Images%22/Sequential_color_palette.png?raw=true)

Figure 2. Seuential color Pallete

![Figure 3](https://github.com/odecampos33/hackbio-cancer-internship/blob/main/Stage%202%22/Images%22/diverging_color_palette.png?raw=true)
Clustering of genes and samples was performed in three ways:

1. Clustering of genes alone (Figure 2).
2. Clustering of both genes and samples (Figure 3.
3. Clustering of samples alone (Figure 4).

 


 

These visualizations identified distinct gene expression patterns, helping to group related genes and patients based on their molecular signatures.

**Differentially Expressed Genes (DEGs)**

DEGs were identified using a t-statistical test, setting a significance level of 5%. The log2 fold change (FC) was calculated using this formula:

<!--[if gte msEquation 12]><m:oMathPara><m:oMath><i
  style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;font-family:
  "Cambria Math","serif";mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:
  "Times New Roman"'><m:r>LogFC</m:r><m:r>=</m:r></span></i><m:sSub><m:sSubPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
    normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>log</m:r></span></i></m:e><m:sub><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>2</m:r></span></i></m:sub></m:sSub><m:d><m:dPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:dPr><m:e><m:f><m:fPr><span
      style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:
      12.0pt;font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
      mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
      mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
      normal'><m:ctrlPr></m:ctrlPr></span></m:fPr><m:num><i style='mso-bidi-font-style:
      normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
      mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>Mean</m:r><m:r>
       </m:r><m:r>of</m:r><m:r> </m:r><m:r>Group</m:r><m:r>1</m:r></span></i></m:num><m:den><i
      style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
      font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
      mso-bidi-font-family:"Times New Roman"'><m:r>Mean</m:r><m:r> </m:r><m:r>of</m:r><m:r>
       </m:r><m:r>Group</m:r><m:r>2</m:r></span></i></m:den></m:f></m:e></m:d></m:oMath></m:oMathPara><![endif]--><!--[if !msEquation]--><!--[if gte vml 1]><v:shapetype id="_x0000_t75"
 coordsize="21600,21600" o:spt="75" o:preferrelative="t" path="m@4@5l@4@11@9@11@9@5xe"
 filled="f" stroked="f">
 <v:stroke joinstyle="miter"/>
 <v:formulas>
  <v:f eqn="if lineDrawn pixelLineWidth 0"/>
  <v:f eqn="sum @0 1 0"/>
  <v:f eqn="sum 0 0 @1"/>
  <v:f eqn="prod @2 1 2"/>
  <v:f eqn="prod @3 21600 pixelWidth"/>
  <v:f eqn="prod @3 21600 pixelHeight"/>
  <v:f eqn="sum @0 0 1"/>
  <v:f eqn="prod @6 1 2"/>
  <v:f eqn="prod @7 21600 pixelWidth"/>
  <v:f eqn="sum @8 21600 0"/>
  <v:f eqn="prod @7 21600 pixelHeight"/>
  <v:f eqn="sum @10 21600 0"/>
 </v:formulas>
 <v:path o:extrusionok="f" gradientshapeok="t" o:connecttype="rect"/>
 <o:lock v:ext="edit" aspectratio="t"/>
</v:shapetype><v:shape id="_x0000_i1025" type="#_x0000_t75" style='width:175.5pt;
 height:30.75pt'>
 <v:imagedata src="file:///C:\Users\hp\AppData\Local\Temp\msohtmlclip1\01\clip_image001.png"
  o:title="" chromakey="white"/>
</v:shape><![endif]--><!--[if !vml]-->![](file:///C:/Users/hp/AppData/Local/Temp/msohtmlclip1/01/clip_image002.png)<!--[endif]--><!--[endif]-->

Where:

<!--[if gte msEquation 12]><m:oMathPara><m:oMath><m:sSub><m:sSubPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
    normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>μ</m:r></span></i></m:e><m:sub><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>1</m:r></span></i></m:sub></m:sSub><i
  style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;font-family:
  "Cambria Math","serif";mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:
  "Times New Roman"'><m:r>= </m:r></span></i><m:f><m:fPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:fPr><m:num><i style='mso-bidi-font-style:
    normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>1</m:r></span></i></m:num><m:den><m:sSub><m:sSubPr><span
      style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:
      12.0pt;font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
      mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
      mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
      normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
      normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
      mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>n</m:r></span></i></m:e><m:sub><i
      style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
      font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
      mso-bidi-font-family:"Times New Roman"'><m:r>1</m:r></span></i></m:sub></m:sSub></m:den></m:f><m:nary><m:naryPr><m:chr
     m:val="∑"/><m:limLoc m:val="undOvr"/><span style='font-size:12.0pt;
    mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-ascii-font-family:"Cambria Math";mso-fareast-font-family:"Times New Roman";
    mso-hansi-font-family:"Cambria Math";mso-bidi-font-family:"Times New Roman";
    font-style:italic;mso-bidi-font-style:normal'><m:ctrlPr></m:ctrlPr></span></m:naryPr><m:sub><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>i</m:r><m:r>=1</m:r></span></i></m:sub><m:sup><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>n</m:r><m:r>1</m:r></span></i></m:sup><m:e><m:sSub><m:sSubPr><span
      style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:
      12.0pt;font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
      mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
      mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
      normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
      normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
      mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>x</m:r></span></i></m:e><m:sub><i
      style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
      font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
      mso-bidi-font-family:"Times New Roman"'><m:r>i</m:r></span></i></m:sub></m:sSub></m:e></m:nary><m:d><m:dPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:dPr><m:e><i style='mso-bidi-font-style:
    normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>Mean</m:r><m:r>
     </m:r><m:r>of</m:r><m:r> </m:r><m:r>Group</m:r><m:r>1</m:r></span></i></m:e></m:d><i
  style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;font-family:
  "Cambria Math","serif";mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:
  "Times New Roman"'><m:r> ; </m:r></span></i><m:sSub><m:sSubPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
    normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>x</m:r></span></i></m:e><m:sub><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>i</m:r></span></i></m:sub></m:sSub><i
  style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;font-family:
  "Cambria Math","serif";mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:
  "Times New Roman"'><m:r>t</m:r><m:r>h</m:r><m:r>e</m:r><m:r> </m:r><m:r>individual</m:r><m:r>
   </m:r><m:r>values</m:r><m:r> </m:r><m:r>in</m:r><m:r> </m:r><m:r>group</m:r><m:r>1</m:r></span></i></m:oMath></m:oMathPara><![endif]--><!--[if !msEquation]--><!--[if gte vml 1]><v:shape id="_x0000_i1025"
 type="#_x0000_t75" style='width:363.75pt;height:41.25pt'>
 <v:imagedata src="file:///C:\Users\hp\AppData\Local\Temp\msohtmlclip1\01\clip_image003.png"
  o:title="" chromakey="white"/>
</v:shape><![endif]--><!--[if !vml]-->![](file:///C:/Users/hp/AppData/Local/Temp/msohtmlclip1/01/clip_image004.png)<!--[endif]--><!--[endif]-->

<!--[if gte msEquation 12]><m:oMathPara><m:oMath><m:sSub><m:sSubPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
    normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>μ</m:r></span></i></m:e><m:sub><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>2</m:r></span></i></m:sub></m:sSub><i
  style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;font-family:
  "Cambria Math","serif";mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:
  "Times New Roman"'><m:r>=</m:r></span></i><m:f><m:fPr><span style='font-size:
    12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;font-family:
    "Cambria Math","serif";mso-ascii-font-family:"Cambria Math";mso-fareast-font-family:
    "Times New Roman";mso-hansi-font-family:"Cambria Math";mso-bidi-font-family:
    "Times New Roman";font-style:italic;mso-bidi-font-style:normal'><m:ctrlPr></m:ctrlPr></span></m:fPr><m:num><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>1</m:r></span></i></m:num><m:den><m:sSub><m:sSubPr><span
      style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:
      12.0pt;font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
      mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
      mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
      normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
      normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
      mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>n</m:r></span></i></m:e><m:sub><i
      style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
      font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
      mso-bidi-font-family:"Times New Roman"'><m:r>2</m:r></span></i></m:sub></m:sSub></m:den></m:f><m:nary><m:naryPr><m:chr
     m:val="∑"/><m:limLoc m:val="undOvr"/><span style='font-size:12.0pt;
    mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-ascii-font-family:"Cambria Math";mso-fareast-font-family:"Times New Roman";
    mso-hansi-font-family:"Cambria Math";mso-bidi-font-family:"Times New Roman";
    font-style:italic;mso-bidi-font-style:normal'><m:ctrlPr></m:ctrlPr></span></m:naryPr><m:sub><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>j</m:r><m:r>=1</m:r></span></i></m:sub><m:sup><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>n</m:r><m:r>2</m:r></span></i></m:sup><m:e><m:sSub><m:sSubPr><span
      style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:
      12.0pt;font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
      mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
      mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
      normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
      normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
      mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>y</m:r></span></i></m:e><m:sub><i
      style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
      font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
      mso-bidi-font-family:"Times New Roman"'><m:r>j</m:r></span></i></m:sub></m:sSub></m:e></m:nary><m:d><m:dPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:dPr><m:e><i style='mso-bidi-font-style:
    normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>Mean</m:r><m:r>
     </m:r><m:r>of</m:r><m:r> </m:r><m:r>Group</m:r><m:r>2</m:r></span></i></m:e></m:d><i
  style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;font-family:
  "Cambria Math","serif";mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:
  "Times New Roman"'><m:r> ; </m:r></span></i><m:sSub><m:sSubPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
    normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>y</m:r></span></i></m:e><m:sub><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>j</m:r><m:r> </m:r></span></i></m:sub></m:sSub><i
  style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;font-family:
  "Cambria Math","serif";mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:
  "Times New Roman"'><m:r>t</m:r><m:r>h</m:r><m:r>e</m:r><m:r> </m:r><m:r>individual</m:r><m:r>
   </m:r><m:r>values</m:r><m:r> </m:r><m:r>in</m:r><m:r> </m:r><m:r>group</m:r><m:r>2</m:r></span></i></m:oMath></m:oMathPara><![endif]--><!--[if !msEquation]--><!--[if gte vml 1]><v:shape id="_x0000_i1025"
 type="#_x0000_t75" style='width:364.5pt;height:42.75pt'>
 <v:imagedata src="file:///C:\Users\hp\AppData\Local\Temp\msohtmlclip1\01\clip_image005.png"
  o:title="" chromakey="white"/>
</v:shape><![endif]--><!--[if !vml]-->![](file:///C:/Users/hp/AppData/Local/Temp/msohtmlclip1/01/clip_image006.png)<!--[endif]--><!--[endif]-->

Then:__<!--[if gte msEquation 12]><m:oMathPara><m:oMath><i
  style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;font-family:
  "Cambria Math","serif";mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:
  "Times New Roman"'>LogFC<m:r>=</m:r></span></i><m:sSub><m:sSubPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i style='mso-bidi-font-style:
    normal'><span style='font-size:12.0pt;font-family:"Cambria Math","serif";
    mso-fareast-font-family:"Times New Roman";mso-bidi-font-family:"Times New Roman"'><m:r>log</m:r></span></i></m:e><m:sub><i
    style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
    mso-bidi-font-family:"Times New Roman"'><m:r>2</m:r></span></i></m:sub></m:sSub><m:d><m:dPr><span
    style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:12.0pt;
    font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
    mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
    mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
    normal'><m:ctrlPr></m:ctrlPr></span></m:dPr><m:e><m:f><m:fPr><span
      style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:
      12.0pt;font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
      mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
      mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
      normal'><m:ctrlPr></m:ctrlPr></span></m:fPr><m:num><m:sSub><m:sSubPr><span
        style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:
        12.0pt;font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
        mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
        mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
        normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i
        style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
        font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
        mso-bidi-font-family:"Times New Roman"'><m:r>μ</m:r></span></i></m:e><m:sub><i
        style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
        font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
        mso-bidi-font-family:"Times New Roman"'><m:r>1</m:r></span></i></m:sub></m:sSub></m:num><m:den><m:sSub><m:sSubPr><span
        style='font-size:12.0pt;mso-ansi-font-size:12.0pt;mso-bidi-font-size:
        12.0pt;font-family:"Cambria Math","serif";mso-ascii-font-family:"Cambria Math";
        mso-fareast-font-family:"Times New Roman";mso-hansi-font-family:"Cambria Math";
        mso-bidi-font-family:"Times New Roman";font-style:italic;mso-bidi-font-style:
        normal'><m:ctrlPr></m:ctrlPr></span></m:sSubPr><m:e><i
        style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
        font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
        mso-bidi-font-family:"Times New Roman"'><m:r>μ</m:r></span></i></m:e><m:sub><i
        style='mso-bidi-font-style:normal'><span style='font-size:12.0pt;
        font-family:"Cambria Math","serif";mso-fareast-font-family:"Times New Roman";
        mso-bidi-font-family:"Times New Roman"'><m:r>2</m:r></span></i></m:sub></m:sSub></m:den></m:f></m:e></m:d></m:oMath></m:oMathPara><![endif]--><!--[if !msEquation]--><!--[if gte vml 1]><v:shape id="_x0000_i1025"
 type="#_x0000_t75" style='width:99pt;height:28.5pt'>
 <v:imagedata src="file:///C:\Users\hp\AppData\Local\Temp\msohtmlclip1\01\clip_image007.png"
  o:title="" chromakey="white"/>
</v:shape><![endif]--><!--[if !vml]-->![](file:///C:/Users/hp/AppData/Local/Temp/msohtmlclip1/01/clip_image008.png)<!--[endif]--><!--[endif]-->__

The cut-offs were >1.5 for upregulation and <-1.5 for downregulation. From the analysis, 10 genes were upregulated and 2 were downregulated (Figure 5). A volcano plot was created to visualize these results (Figure 5).

**Functional Enrichment Analysis**

Using ShinyGO (Ge et al., 2019b), functional enrichment analysis was performed to identify key biological pathways affected by the upregulated and downregulated genes (Figures 6 and 7). Lollipop plots were created to show the top five pathways, scaling the points according to the negative log10 of the p-value to reflect their significance.

The top 3 upregulated pathways include:

1. **Glutathione derivative metabolic process**: Critical in neutralizing reactive oxygen species (ROS) and promoting cancer cell survival through redox balance (Kennedy et al., 2020).
2. **Glutathione derivative biosynthetic process**: Involves increased glutathione production, allowing glioblastoma cells to survive oxidative stress (Backos et al., 2012).
3. **Proteolysis**: Supports extracellular matrix degradation and enhances tumor growth by remodeling the tumor microenvironment (Bischof et al., 2017).

**Top 3 Downregulated Pathways:**

1. **Ribosomal small subunit assembly**: Ribosomal proteins support tumor progression by maintaining stem-cell-like properties in glioblastoma (Hide et al., 2022).
2. **Maturation of SSU-rRNA from tricistronic rRNA transcript**: Disruptions in ribosomal function promote carcinogenesis (Shirakawa et al., 2021).
3. **Maturation of SSU-rRNA**: Linked to ribosomal disorders, contributing to tumor development (McElreavey et al., 2022).

These analyses provided insights into the pathways driving glioblastoma pathogenesis, offering potential targets for therapeutic intervention.

 


# **REFERENCES**

Backos, D.S., Franklin, C.C. and Reigan, P., 2012. The role of glutathione in brain tumor drug resistance. _Biochemical Pharmacology_, 83(8), pp.1005-1012. doi: 10.1016/j.bcp.2011.11.016.

Bischof, J., Westhoff, M.-A., Wagner, J.E., et al., 2017. Cancer stem cells: The potential role of autophagy, proteolysis, and cathepsins in glioblastoma stem cells. _Tumor Biology_, 39(3). doi:10.1177/1010428317692227.

Brosseau, S., Viala, M., Varga, A., Planchard, D., Besse, B. and Soria, J.C., 2015. 3rd generation's TKI in lung cancer non-small cell EGFR-mutated having acquired a secondary T790M resistance. _Bulletin du Cancer_, 102, pp.749–757. doi: 10.1016/j.bulcan.2015.05.001.

Brown, N.F., Ottaviani, D., Tazare, J., Gregson, J., Kitchen, N., Brandner, S., Fersht, N. and Mulholland, P., 2022. Survival outcomes and prognostic factors in glioblastoma. _Cancers_, 14(13). doi: 10.3390/cancers14133161.

Chen, J., Zeng, F., Forrester, S.J., Eguchi, S., Zhang, Z. and Harris, R.C., 2016. Expression and function of the epidermal growth factor receptor in physiology and disease. _Physiological Reviews_. doi: 10.1152/PRV-00030-2015.

Gilard, V., Tebani, A., Dabaj, I., Laquerrière, A., Fontanilles, M., Derrey, S., Marret, S. and Bekri, S., 2021. Diagnosis and management of glioblastoma: A comprehensive perspective. _Journal of Personalized Medicine_, 11(4). doi: 10.3390/jpm11040258.

Hide, T., Shibahara, I., Inukai, M., Shigeeda, R. and Kumabe, T., 2022. Ribosomes and ribosomal proteins promote plasticity and stemness induction in glioma cells via reprogramming. _Cells_, 11(14), p.2142.

Huang, W., Li, J., Zhu, H., Qin, X., Chen, C., Wang, B., Wei, J., Song, Y., Lu, X., Li, Z., et al., 2023. A novel EGFR variant EGFRx maintains glioblastoma stem cells through STAT5. _Neuro-Oncology_, 26, pp.85–99.

Kennedy, L., Sandhu, J.K., Harper, M.E. and Cuperlovic-Culf, M., 2020. Role of glutathione in cancer: From mechanisms to therapies. _Biomolecules_, 10(10), p.1429. doi: 10.3390/biom10101429.

McElreavey, K., Pailhoux, E. and Bashamboo, A., 2022. DHX37 and 46,XY DSD: A new ribosomopathy? _Sexual Development_, 16, pp.194–206.

Singh, S., Barik, D., Lawrie, K., Mohapatra, I., Prasad, S., Naqvi, A.R., Singh, A. and Singh, G., 2023. Unveiling novel avenues in mTOR-targeted therapeutics: Advancements in glioblastoma treatment. _International Journal of Molecular Sciences_, 24, p.14960.

Shirakawa, Y., Ohta, K., Miyake, S., Kanemaru, A., Kuwano, A., Yonemaru, K., Uchino, S., Yamaoka, M., Ito, Y., Ito, N., Hide, T., Shinojima, N., Mukasa, A., Saito, H. and Jono, H., 2021. Glioma cells acquire stem-like characters by extrinsic ribosome stimuli. _Cells_, 10(11), p.2970. doi: 10.3390/cells10112970.

Ge, S. X., Jung, D., & Yao, R. (2020). ShinyGO: a graphical gene-set enrichment tool for animals and plants. _Bioinformatics (Oxford, England)_, _36_(8), 2628–2629. doi:10.1093/bioinformatics/btz931.
