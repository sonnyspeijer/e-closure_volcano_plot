# e-Closure Volcano Plot

This application implements an interactive volcano plot with simultaneous false discovery rate (FDR) control, based on methods constructed using the e-closure principle. The online version is accessible [here](https://sonnyspeijer.shinyapps.io/e-closure_volcano_plot/). 

## Motivation

A common practice in genomics is to adjust p-values using the Benjamini-Hochberg (BH) procedure to control the FDR when testing multiple hypotheses at the same time. If this procedure yields a large discovery set, a researcher may want to reduce its size. This can be accomplished using the volcano plot, which filters the results not only by p-value but also by effect size, usually expressed as a fold change.

The issue with this approach is that applying a fold change filter to a discovery set with FDR control does not guarantee that the resulting subset also has FDR control. This can be circumvented by constructing an FDR method using the e-closure principle, which returns a collection of discovery sets with <i>simultaneous</i> FDR control. Simultaneous control means that the filters can be changed post hoc without invalidating the results. Consequently, researchers can retain FDR control as long as the filtered discovery set is contained within the collection returned by the e-closure principle.
