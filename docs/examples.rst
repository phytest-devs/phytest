==============
Examples
==============

The Phytest organisation of GitHub contains serval example repositories (https://github.com/phytest-devs) that show how Phytest can integrate into standard phylogenetic analyses and scenarios.

Nextstrain Example
------------------

The Nextstrain pipeline is widely used for pathogen phylogenetic analysis.
In this example we modify the Nextstrain zika-tutorial (https://github.com/nextstrain/zika-tutorial)
to include testing with Phytest. This modified pipeline is available as an example repository in the phytest-devs GitHub organisation
https://github.com/phytest-devs/phytest-nextstrain-example and provides an example of using Phytest for quality control in a Snakemake
pipeline.

Phytest is included in the pipeline to ensure the alignment and maximum likelihood tree meet explicit quality
requirements before proceeding though the pipeline. Only if all the tests pass will the pipeline continue, thus savings computational resources.
The resulting HTML report provides details of any failed tests so that the offending data can be removed.
While Augur (the Nextstrain toolkit) has some ability to refine/filter tree and alignment files,
Phytest adds highly a customizable testing framework to the pipeline that ensures the quality of the analysis.

Temporal Signal Example
-----------------------

A repository containing the code for this example can be found at https://github.com/phytest-devs/phytest-temporal-signal-example.

Temporal signal in an important prerequisite to many Bayesian phylogenetic analyses. In this example we use Phytest to ensure the
data-set meets the minimum temporal signal requirements for Bayesian analyses. Temporal signal analysis can help to detect
problematic sequences and potential issues before heading on to a Bayesian phylogenetic analysis e.g. with BEAST.
Here, we use data from from the TempEst tutorial https://beast.community/tempest\_tutorial.

TempEst is a useful program for performing temporal signal analysis, however, it is not possible to easily automate the TempEst graphical user interface.
Internally, Phytest uses TimeTree to perform a root-to-tip regression, allowing users to automate temporal signal testing.
The :code:`Tree.assert_root_to_tip` method is used for testing temporal signal and provides arguments for testing the
coefficient of determination, estimated rate and root date. The Phytest Tree class also implements methods for exploring and plotting results.

Continuous Testing Example
--------------------------

In this example Phytest is used to test data shared on GitHub every time the data is updated (https://github.com/phytest-devs/phytest-continuous-testing-example).

Tests are run against the phylogenetic data using the Continuous Integration features that are freely available
through GitHub (other services are also available). Using Phytest through GitHub Actions (https://github.com/features/actions)
ensures that anytime the data changes (common during development) they still meet the requirements defined in the tests.
