# Final Project for BF591 (R for Biological Sciences)

## Project Overview
The RShiny app was developed as part of the final project and the aim was to allow exploration a differential expression dataset comparing post-mortem Huntington’s Disease Brain data. The data for this is publically available through [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810). The paper used for reference is also [publically available](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4670106/). The theme and design for the app is entirely done using CSS and not by utilising in-built RShiny theme libraries.

The different tabs in the app correspond to the data used as input. An explanation of the files used which are present in the data directory:

- GSE64810_series_matrix.txt - Metadata file which contains various information about the samples/experiments.
- GSE64810_mlhd_DESeq2_norm_counts_adjust.csv - Normalised counts adjusted DESeq data.
- GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.csv and txt - Differentially Expressed with outliers removed data.
- human_id2gene.txt - Used for conversion between Ensembl IDs.
- fgsea_results.csv - Results of fgsea analysis which was performed using original script fgsea_analysis.R present in the main directory.
- c3.tft.v2023.2.Hs.symbols.gmt - Hallmark gene sets which was used in fgsea analysis.

## Tabs in the app

### 1. Sample Metadata
   
This tab is used to display the sample metadata information. It uses the series matrix txt file as input and all the columns are filterable and there is a search bar which can be used. There are 3 additional nested tabs:
* Summary: This tab shows the column information as well as the data type of the column. It also takes mean of the column if the data type is numerical in nature.
* Metadata: This displays the entire metadata information present in the input txt. Horizontal scrolling is enabled for ease of use.
* Violin Plot: This tab dynamically displays violin plots between continous variables of the users choosing. The options between the columns that can be plotted is given by usind a drop down.

Metadata tab:
![Screenshot (31)](https://github.com/posaumya/RShiny_Project/assets/144373823/0f01299a-8072-47a5-845a-aac23020a95d)


### 2. Counts Matrix
This tab is used to display the counts matrix information which is obtained by using the normalised counts csv file. There are 4 nested tabs which use information returned from 2 sliders which are used to dynamically take user input in gene percentile of variance and for choosing genes that have samples with non-zero values. The additional nested tabs are:
* Filtered Data: This tab shows a summary of genes/samples passing/failing the user input threshold of cutoff (slider values).
* Scatter Plot: This tab displays 2 scatter plots- one is Ranked Median vs Log Variance and the other is Ranked Median vs Number of Zeros. The colours of the data points are based on pass/fail criteria. 
* Heatmap: The heatmap tab shows the heatmap of all samples with their gene expression levels plotted.
* PCA: The PCA tab allows the user to view different number of PCs after PCA is performed in the backend. The input for number of PCs which can be viewed is taken by using a slider input. The percentage contribution of each PC is dynamically displayed in the legend.

Scatter Plot:
![Screenshot (30)](https://github.com/posaumya/RShiny_Project/assets/144373823/f60c0912-23ac-4fa0-bafa-234d95d3ac7a)


### 3. Differential Expression
This tab uses the differentially expression genes csv as input and has 2 additional nested tabs:
* Metadata: This tab displays the entire csv metadata information. Horizontal scrolling is enabled for ease of use, and all columns are filterable with a search bar present for faster pulling of information.
* Plot: This tab has 2 additional nested tabs which are volcano plot and filtered table. The user chooses a padj threshold value as a cutoff by selecting a value using a slider, and can also change the variables which are plotted on the x and y axis by choosing column names via radio buttons. The data points which pass the padj filtering are plotted in different colours, which can also be dynamically changed by users. The data points which pass this filter are present in the filtered table tab.

Volcano Plot:
![image](https://github.com/posaumya/RShiny_Project/assets/144373823/51e46263-68be-4cdf-8b6b-9cb5e56d53fa)

### 4. GSEA
The last tab is used to explore the fgsea analysis results which was performed using the R script fgsea_analysis.R present in the main directory. There are 3 nested tabs:
* Barplot: This tab plots a barplot of the Human MSig C3 TFT Gene Sets using Normalised Enrichment Scores (NES) on the X-axis. Pathways are picked based on user input of padj value which is obtained using a slider.
* Table: This tab also takes padj threshold input from the user and filters pathways based on the value. The pathways can be further filtered based on if they are positive/negative or all pathways can be displayed using radio buttons. This filtered table can be downloaded using a download button for further analysis. Again horizontal scrolling is enabled for ease of use, columns are filterable and there is search bar for fast retrieval.
* Scatter Plot: This tab has a scatter plot of NES vs -log10(padj) values. The padj values are taken as input from the user using a slider, and the values which pass are highlighted.

Barplot:
![image](https://github.com/posaumya/RShiny_Project/assets/144373823/f0595a7a-6b6e-4a46-bdaa-419430c8277e)

Scatter Plot:
![image](https://github.com/posaumya/RShiny_Project/assets/144373823/aaf932f4-8752-4f25-9629-619527cdefb0)





