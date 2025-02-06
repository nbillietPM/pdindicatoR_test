library(ape)
library(ggplot2)

angiospermTree2012 <- read.tree("./phylogeneticAngiosperm/DaPhnE_01/DaPhnE_01.tre")
angiospermTree2012$tip.label <- lapply(angiospermTree2012$tip.label, function(label) {gsub("_", " ", label)})
"
Phylogenetic tree with 5122 tips and 3345 internal nodes.

Tip labels:
  Huperzia_selago, Lycopodiella_inundata, Lycopodium_clavatum, Lycopodium_annotinum, Diphasiastrum_zeilleri, Diphasiastrum_tristachyum, ...
Node labels:
  Vascular_plants, Lycopodiophyta, Lycopodiales, N4, N6, Lycopodium, ...

Rooted; includes branch length(s).
"
#The newer tree is in a Bayesian format and requires to be read as a nexus format
angiospermTree2020 <- read.nexus("./phylogeneticAngiosperm/oo_330891/oo_330891.tre")
angiospermTree2020$tip.label <- lapply(angiospermTree2020$tip.label, function(label) {gsub("_", " ", label)})
"
Phylogenetic tree with 36106 tips and 36104 internal nodes.

Tip labels:
  Amborella_trichopoda, Trithuria_bibracteata, Trithuria_submersa, Trithuria_austinensis, Trithuria_australis, Trithuria_filamentosa, ...

Unrooted; includes branch length(s).
"

