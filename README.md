
Documentation at [splicing.cengen.org](http://splicing.cengen.org).

# DAS by neuron

Shiny app to help make comparisons of alternative splicing between neurons on the CeNGEN data.

The local quantification was performed with MAJIQ, and its results can be examined using the online VOILA app at [http://splicingapps.cengen.org/voila](http://splicingapps.cengen.org/voila). However, this provides one page per gene, and does not provide functionality to discover genes varying between neurons.

This app allows the user to enter neurons of interest and gets events and genes that appear differentially alternativelty spliced (DAS).

## Available usage

In the "neuron pair" tab, one can input two neurons, and we use the statistics from MAJIQ to select events that are DAS.

In the "neuron sets" tab, one can input two sets of neurons. We then recover the PSI of each junction in neurons of these sets and perform a t-test. As the t-test tends to be robust, the results will usually point to set of DAS events for further analysis. However, depending on the exact neuron sets used as input, the test may fail to correctly control false discoveries. Thus, the resulting list of DAS event should not be considered the "true" set of DAS events and may be misleading. In particular, any conclusion such as "this set of neurons differs more from these neurons than those neurons" is not valid.

## Event names

MAJIQ provides event IDs of the form "WBGene00006064:t:5311643-5311763". These identifiers are unique and robust, and should be the main way to describe an event. However, they are unpractical to reason with and cross-reference between several neuron types. For convenience, we also attribute a name to each event, e.g. "Quantasha". These names were randomly attributed to each event, and may change in future updates, but they are consistent within the application. They are names that were given to at least 5 babies in the United States between 1880 and 2017, as provided by the package [`{babynames}`](https://hadley.github.io/babynames/).

