## Tutorial 4

**Objective:** Perform a pathway enrichment analysis using your DESeq2 results.

The next step in our analyses will be to functionally annotate our genes, assign them to higher-level categories, and determine whether any categories have an over-representation of differentially expressed genes.
[Prokka](https://github.com/tseemann/prokka) and [bakta](https://github.com/oschwengers/bakta) are commonly used to annotate bacterial genomes (note that prokka is no longer maintained). Higher-level classification can be achieved with various tools, such as the Gene Ontology (GO) database [https://geneontology.org/](https://geneontology.org/). Genes are hierarchically classified into three main sub-ontologies: Molecular function, Biological process, or Cellular component as well as more specialized related terms.

[KEGG](https://www.kegg.jp/) is another resource for assigning annotation and predicting higher-level functions.  Genes or proteins are annotated with KO (KEGGG orthology) identifiers or K numbers, which are mapped to molecular networks.  KOs are characterized via KEGG Pathway and Brite hierarchies.  KEGG pathways are hierarchically organized into the following main categories:
* Metabolism
* Genetic Information Processing
* Environmental Information Processing
* Cellular Processes
* Organismal systems
* Human Diseases
* Drug Development

These main categories are further subdivided in sub-categories - such as "Carbohydrate" and "Energy" Metabolism -  which are then broken down in to pathways, such as "Citrate acid (TCA) cycle" and "Oxidative phosphorylation."

We will use KEGG to annotate our genes (assign K numbers) and [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler) to perform enrichment analyses based on KEGG classifications.

Enrichment analyses identify pathways or gene categories that have a higher proportion of differentially expressed genes than by chance alone.
We calculate pathway enrichment by dividing the proportion of DEGs in our pathway (DEGs in pathway/total number of genes in pathway) by the proportion of DEGs in our dataset (total number of annotated DEGs/total number of annotated genes) or our background gene set.  
We can test the null hypothesis with the hypergeometric test that determines the probability of getting k or more number of DEGs inside a gene set based on the hypergeometric distribution (a discrete probability distribution that describes the number of successes in a series
of draws from a finite population without replacement). **Note that the hypergeometric test is equivalent to a one-tailed Fisher's exact test.**

<br>

# Generate a vector of differentially expressed genes based on your DESeq2 results

Summarize the results from the DEG analysis

	summary(res, alpha=0.05)

> This quickly allows us to see how many genes are up or down-regulated at an alpha value of 0.05.

<br>

Convert the res object to a dataframe for easier manipulation.

	res_table = as.data.frame(res)

<br>

Filter the table based on adjusted p-value and log2fold change using base R functions.

	sigDE = res_table[which(res_table$padj <= 0.05 & abs(res_table$log2FoldChange) >=1),]
	
<br>

We can also filter with functions from the `dplyr` package (part of the [tidyverse](https://tidyverse.org/) package installation).

	sigDE <- res_table %>%
  		filter(padj <= 0.05, abs(log2FoldChange) >= 1)
  		
<br>

We need a list of genes (really what R calls a vector) for pathway enrichment analysis with `clusterProfiler` package.

Conveniently we can generate this from the row names of our `sigDE` table.

	genes = rownames(sigDE)
	is.vector(genes)
	
> **`is.vector()`** checks the format of our delist object.  It is indeed a vector.

<br>

# Search for the appropriate genome and its annotation in the KEGG database.

To run a KEGG pathway enrichment analysis in `clusterProfiler`, we need to associate the locus tags with KEGG's functional annotation.

KEGG genomes are referenced with three character codes.  We can search for those associated with _Pseudomonas aeruginosa_.  

**The search is case sensitive.**

	search_kegg_organism("Pseudomonas aeruginosa")
	
<br>

# Perform an enrichment analysis with your vector of DE genes.

Now we can select the appropriate code, which is `pau` for `_Pseudomonas aeruginosa_ UCBPP-PA14` and perform our enrichment analysis.

	en = enrichKEGG(gene=genes, organism='pau')

The results are in a dataframe with columns for:
* category: names of enriched KEGG categories
* subcateogry: names of enriched KEGG subcategories
* Description: names of enriched KEGG pathways
* GeneRatio: number of DE genes in pathway to total number of genes in pathway
* BgRatio: background ration of the number of total DE genes to the total number of genes in dataset
* FoldEnrichment: ratio of GeneRation to BgRatio
* pvalue: p-value
* p.adjust: BH adjusted p-value

<br>

# Visualize the results with functions from [enrichplot](https://github.com/YuLab-SMU/enrichplot) package.

Note that `enrichplot` is loaded with [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler).

Visualize the results as a barplot

	barplot(en)

> The DE counts for the first 10 enriched categories are shown along the x-axis. <br>
> Bars are colored by the strength of the adjusted p-value. <br>

<br>

Visualize the results as a bubbleplot.

	dotplot(en)
	
> Here, the gene ratio (DEGs in pathway/annotated genes in pathway) is represented along the x-axis and the size of the dot corresponds to the number of DE genes in an enriched KEGG category. <br>
> Again, the color indicates the strength of the adjusted p-value. <br>
	
* Are these results consistent with our expectations given the role of pqsE in regulating quorum sensing and biofilm formation? <br>
* The default KEGG IDs are the locus tags for bacteria.  We could convert our locus tags to KOs and see if our enrichment results differ. <br>

<br>

We can convert our locus tags into different gene identifiers with the `bitr_kegg()` function. <br>

	bitr_kegg(genes, fromType="kegg", toType="uniprot", organism="pau", drop = TRUE)
	
> Unfortunately, this tool does not allow us to convert between locus tags and KO identifiers (K numbers) of KEGG.

<br>
	
To convert our locus tags to K numbers, we need to download the KEGG file for our organism of interest from the [KEGG website](https://www.kegg.jp/kegg/).

* Navigate to KEGG website and scroll down to the `KEGG organisms` link under the `KEGG database` options.  Since we know the KEGG organism code (`pau`), we can type it into the box and click `go`. 
* This brings us to the genome page for `Pseudomonas aeruginosa UCBPP-PA14`.  Click on the link for `Brite hierarchy` at the top of the page. 
* This brings us to a page for the BRITE functional hierarchies, including a link to the `KEGG orthology` for our organism of interest.  Clicking this link will bring you to a page for the KEGG Orthology (KO) for Pseudomonas aeruginosa UCBPP-PA14. 
* Click the right-most down arrow to expand the contents of the page, which shows you the locus tag mapped to its gene abbreviation, gene annotation, and associated K number.  
* You can download this page to a text file by selecting the `Download htext`.  With a little reformatting, you can use this file to convert locus tags to K identifiers.

<br>

For now, this step has been done for you.

	ko_list = c("K00077","K00370","K00371","K00373","K00374","K00376","K00763","K01028","K01141","K01183","K01626","K01738","K01895","K01980","K02047","K02305","K02411","K02415","K02575","K02651","K03638","K04561","K06998","K08738","K11473","K13063","K15864","K19341","K20258","K20260","K20261","K20262","K21103","K22225","K24843","K24844","K24867")
	
	en_ko = enrichKEGG(gene=ko_list, organism="ko")
	dotplot(en_ko)

* Did our results change?

<br>

# Pull a pathway map from KEGG and convert to a graph network with the [KEGGgraph](https://github.com/Accio/KEGGgraph) package.

	tmp = "pau02024.xml"
	retrieveKGML(pathwayid="02024", organism='pau', destfile='tmp')
	g = parseKGML2Graph("tmp", expandGenes=TRUE)

> We retrieve the KGML file for Pseudomonas quorum sensing and parse the information to generate a graph object. <br>
> We can plot the KEGGgraph graph or we can convert it to an `igraph` graph for additional manipulation. <br>

<br>

# Convert the graph to an [igraph](https://igraph.org/) graph object.

	igraph_object = graph_from_graphnel(g)
	plot(igraph_object, layout=layout_nicely(igraph_object))
	
> The figure represents the quorum sensing pathway as a series of nodes or vertices connected by edges. <br>
> The vertices represent the proteins that comprise the pathway, the edges represent the interactions among proteins. <br>
> In this case, our graph is directed, with the arrows showing the direction of the interaction.

<br>

# Resize vertices based on their connectedness. 

The connectedness among nodes/vertices represents the number of edges connecting them to other vertices.

	V(igraph_object)$degree = degree(igraph_object)
	plot(igraph_object, layout=layout_nicely(igraph_object), vertex.size=V(igraph_object)$degree, vertex.label=NA)
	
> To make this graph more informative, we can resize and color the graph's attributes (i.e. nodes and edges) by data from our DESeq2 and enrichment analyses.

<br>

# Color vertices based on their differential expression.

If you recall, we extracted the names of our differentially expressed genes from the DESeq2 results table, `res_table` and saved them to the `genes` vector.
We can color the corresponding vertices in our igraph object using this information.
First ensure that our gene names from our DESeq2 results match the vertex names. In this case, the vertex names are displayed as the locus tag, preceeded by `pau:`. 

We can change the names of our vertices so that they match the gene names of our DESeq2 results by dropping the `pau:` from the beginning of each locus tag:

	V(igraph_object)$name = gsub("pau:", "", V(igraph_object)$name)
	
> **`gsub`** is a search and replace function that accepts regular expressions. <br>
> We overwrite the vertex names of our `igraph_object` object, searching for instances of `pau:` and replacing the expression with nothing or `""`.
<br>

Now we can use an `if else` statement to color any vertex names present in our `genes` list.
	
	V(igraph_object)$color = ifelse(V(igraph_object)$name %in% genes, "lightblue", "gray")
	plot(igraph_object, layout=layout_nicely(igraph_object), vertex.size=V(igraph_object)$degree, vertex.label=NA, vertex.color=V(igraph_object)$color)
	
> We populate a `color` attribute for our graph vertices depending on whether the names of the vertices are found in `genes`. <br>
> The `ifelse` statement has the following syntax: if the vertex names of `igraph_object` are in `genes`, then assign `lightblue`, else assign `gray`.

<br>

# Label only vertices that are differentially expressed.

We can use another `ifelse` statement to label only certain nodes on our network graph and minimize clutter.

	plot(igraph_object, layout=layout_as_tree(igraph_object), vertex.color=V(igraph_object)$color, vertex.label=ifelse(V(igraph_object)$name %in% genes, V(igraph_object)$name, NA))

<br>

# Resize vertices based on their log2FoldChanges.

Create a dataframe with the vertex names.

	V_names <- as.data.frame(V(igraph_object)$name)
	colnames(V_names) <- "name"

<br>

Add a corresponding column, named `name` in our DESeq2 results table that contains the locus tag information.

	res_table$name <- rownames(res_table)

<br>

We can now join the information of the two tables together based on the matching locus tags in both "name" columns.
Only the data from the locus tags that match in `V_names` and `res_table` are saved to the `vertex_attributes` dataframe.

	vertex_attributes = left_join(V_names, res_table, by = "name")

<br>

Now store the absolute log2FoldChange values as the size attribute of our graph vertices.

	V(igraph_object)$size = abs(vertex_attributes$log2FoldChange)
	plot(igraph_object, layout=layout_as_tree(igraph_object), vertex.color=V(igraph_object)$color, vertex.label=ifelse(V(igraph_object)$name %in% genes, V(igraph_object)$name, NA), vertex.size=V(igraph_object)$size)

<br>

# Color vertices based on up or down-regulation of gene expression.

Make a vertex attribute for log2FoldChanges and color the nodes depending on its sign.

	V(igraph_object)$lfc = vertex_attributes$log2FoldChange
	V(igraph_object)$color = ifelse(V(igraph_object)$lfc > 0, "red", "blue")
	plot(igraph_object, vertex.color=V(igraph_object)$color, vertex.label=ifelse(V(igraph_object)$name %in% genes, V(igraph_object)$name, NA), vertex.size=V(igraph_object)$size, vertex.label.dist=-1, vertex.label.cex=0.5, vertex.label.degree=pi/2, edge.arrow.size=0.5)

> All nodes are colored regardless of whether they are up or down-regulated, although only nodes with an adjusted p-value < 0.05 are labeled.

<br>

# Conditionally color vertices blue, red, or gray depending on significance and value of change.

Make a vector `vertex_colors` with empty entries (the number of entries corresponding to the number of vertices in our graph), then populate it with different color assignments based on filtering the `vertex_attributes` dataframe.

	vertex_colors <- character(nrow(vertex_attributes))
	vertex_colors[vertex_attributes$padj <= 0.05 & vertex_attributes$log2FoldChange <= -1] = "blue"
	vertex_colors[vertex_attributes$padj <= 0.05 & vertex_attributes$log2FoldChange >= 1] = "red" # there happen to be no genes that fit these criteria
	vertex_colors[vertex_attributes$padj > 0.05 | abs(vertex_attributes$log2FoldChange) < 1] = "gray"

> The first conditional statement searches for entries in the `vertex_attributes` dataframe that have adjusted p-values <= 0.05 *AND* a log2FoldChange <= -1 and assigns the value `blue` to the corresponding row indices of those entries in `vertex_colors`. <br>
> The syntax of second conditional statement mirrors the first but slots the `red` value into `vertex_colors` based on the row indices of entries in `vertex_attributes` that satisfied the specified conditions. <br>
> The third conditional statement identifies entries in vertex_attributes that have an adusted p-value > 0.05 *OR* (`|`) an absolute fold change < 1 (non-differentially expressed genes). <br>

<br>

Genes with adjusted p-values of "NA" do not have a color assigned are represented by empty entries in our color vector.  We can change this with the `gsub` command so that these genes are also colored gray.

	vertex_colors = gsub("^$", "gray", vertex_colors)
	V(igraph_object)$color = vertex_colors
	plot(igraph_object, vertex.color=V(igraph_object)$color, vertex.label=ifelse(V(igraph_object)$name %in% genes, V(igraph_object)$name, NA), vertex.size=V(igraph_object)$size, vertex.label.dist=-1, vertex.label.cex=0.5, vertex.label.degree=pi/2, edge.arrow.size=0.5)

<br>

You can also create a sparser matrix by deleting nodes or edges.

	g = delete_edges(igraph_object, E(mc_graph)[weight<4])
	g = delete_vertices(igraph_object, V(mc_graph)[degree<1])

