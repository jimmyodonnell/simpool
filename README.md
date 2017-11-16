### Is environmental DNA reliable?

Environmental DNA (eDNA) is currently used to estimate properties of ecological communities (e.g. richness, evenness, presence of a single species).
However, the accuracy and reliability of those estimate is unknown. 
In contrast to traditional sampling methods (e.g. visual transects), the path from a true ecosystem property to eDNA data consists of a long chain of independent steps. 
(In essence, this amounts to a long sequence of subsampling events. )
Along this chain, there are two fundamental and distinct types of processes that determine the degree to which eDNA data can reliably estimate ecosystem properties: Those that influence the amount of DNA that ends up in an environmental sample, and those that influence the estimation of the type and amount of DNA in that sample. 

**Before environmental sampling: (BES)**
- An organism sheds material that contains DNA.
- That material ends up away from the organism by movement of the organism, the environmental medium, or both.
- DNA degrades over time.


  - **Factors influencing BES:**
    - Biomass and number of individuals and their distance from sample location
    - Shedding: biomass, metabolism, rate of mucous/waste generation, dermal permeability?
    - Movement: speed+direction of organism and environment
    - Degradation: temperature, radiation, enzymatic activity, chemical properties


**After environmental sampling: (AES)**
- DNA degrades after sample collection.
- DNA is isolated from a sample.
- A subset of that DNA is amplified for analysis.
- The relative or absolute amount of target DNA is inferred.


  - **Factors influencing AES:**
    - Time and environment of sample storage after collection
    - Extraneous DNA entering sample after collection
    - Efficiency of sample concentration (e.g. filtration)
    - Time and environment of sample storage after concentration
    - Efficiency of DNA isolation from sample
    - Time and environment of sample storage after concentration (e.g. preservative buffer, number of freeze-thaw cycles)
    - Accuracy of pipetting in each laboratory procedure
    - Measurement accuracy of quantification (e.g. to create standards for qPCR or to pool samples for multiplex sequencing)
    - Efficiency of amplification (template-specific)
    - Measurement accuracy of final sample:
      - qPCR accuracy
      - Sequencing accuracy (e.g. sequencing errors, fragment size bias)

This analysis uses simulation to assess the impact of factors influencing the inference made **after** a sample is collected from the environment. 
Specifically, we focus on:
  - the effect of primer mismatch on sample inference
  - the effect of replicating laboratory procedures on sample inference

#### Primer Mismatch
Mismatch is assumed to be quantifiable over {0,1}, where 0 is complete mismatch and thus no amplification, and 1 is perfect match and thus amplification is not limited by primer mismatch. 
We thus model the distribution of mismatches in a sample as draws from a beta distribution. 
These correspond to situations in which the _variation_ in mismatch among templates is high, medium, low, and absent. 
Each scenario and the parameters of the corresponding beta distribution is represented in the table below.


| Mismatch Variation | Shape1 | Shape2 |
|---|:-:|:-:|
| high   | 1  | 1 |
| medium | 2  | 1 |
| low    | 10 | 1 |
| absent | 1  | 0 |
| real1  | A1  | B1 |
| real2  | A2  | B2 |

Additionally, we quantified the real distribution of mismatches between a set of taxa and each of 2 primer sets. 
((Alternatively, we could just sample from the empirical values themselves.))

We generated data under these conditions:
- True DNA template communities:
  - Community consists of 10 template types (aka species aka sequence variants)
  - Variation in template abundance: all equal, varying, etc. 
  - vary total quantity?
- Variation in in-silico PCR conditions:
  - Variation in template mismatch
  - Variations correspond to true abundance?

Thus, the input to the simulations looks like:

| variable | type | values |
|---|:-:|--:|
| community.id | integer | 1:(Nscenarios \* N_rep?) |
| div.scenario | factor | (even, exp, etc.) |
| species | integer | 1:10 |
| copies.per.ul | integer | 0:? |

The output from the simulations looks like this:

| variable | type | values |
|---|:-:|--:|
| pcr.id | integer | 1:(N_comm \* N_rep \* N_mm ) |
| community.id | integer | 1:N_comm |
| replicate | integer | 1:10 |
| mmv.scenario | factor | hi, med, low, abs, real1, real2 |
| species | integer | 1:10 |
| copies.per.ul | integer | 0:? |

Beyond the bias introduced by PCR mismatch, we assume that measurement of the quantity during or following PCR is perfect (without error).

From this, we can calculate relative abundance in input, and relative abundance in output.
We calculated the error introduced to each template as (OBS - EXP)/EXP.

Finally, because the goal of community eDNA studies is to quantify some property of the community, we calculated richness and evenness, and compared these estimates. 
Because we intentionally use a community with low richness, we focus especially on measures of evenness.
We could convert this first to (OBS - EXP)/EXP in order to have a consistent expectation that points should fall at 0.

```r
center <- function(OBS, EXP){return((OBS - EXP)/EXP)}
boxplot( center(evenness) ~ true_evenness, box.col = levels(mmv.scenario))
```

Finally, for scenarios in which the variation in inferred evenness is high (e.g. mmv = high), we conducted the following analysis:
  - pick two samples from pcr.id
  - calculate difference in evenness -> true
  - get corresponding community.id
  - calculate difference in evenness -> inferred
  - `plot(inferred ~ true, pch = mmv.scenario)`
  - Maybe fit linear model, but perhaps test for heteroscedasticity

#### Replication
Up to this point, we have only considered single replicates of the process. 
However, replication could help overcome biased measurements. 
(Can it? It should only help when bias is introduced by stochastic processes, and here bias is due to fixed primer bias.)


