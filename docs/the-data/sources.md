# Source catalogue

The catalogue lists the source papers, compilations and databases represented in
optmat. Each page reproduces the recorded citation, technique and curation notes,
followed by the dataset identifiers, locations and use designations, recorded
formulas and corrections. It is generated from the committed source records at
every documentation build.

```{toctree}
:maxdepth: 1

../generated/sources/index
```

## What a catalogue entry tells you

The structured entries show which materials and quantities were recorded. The
notes explain contributions that are not captured by those fields: definitions,
interpretations, model assumptions and limitations, for example. A paper may have
formulas but no datasets, or notes without either. A missing DOI means that none
is recorded; it does not imply that the source has no DOI.

The citation is shown as recorded, including abbreviations or incomplete details.
The catalogue does not independently validate bibliographic metadata. A recorded
bibliographic identity status concerns the identification of the paper, not the
scientific correctness of its results.

Dataset metadata and formulas are displayed as structured text. Field meanings,
units and uncertainty conventions are in {doc}`/data` and
{doc}`/the-data/quantities`. The link to the source record provides the complete
TOML, including inline points or the paths of extracted data files.

## From a paper to a model

The designation `primary` makes a dataset eligible for selection; it does not
establish that the dataset enters a recommended model. The fit configuration and
the {doc}`review register </register>` record selection decisions, corrections,
checks and doubts.
Use the relevant scientific chapter in the {doc}`/user-guide/index` to understand
the model's validity and uncertainty.

For a Python example retrieving a dataset and its citation, return to
{doc}`/quickstart`. For the relationship between sources, datasets and models,
read {doc}`/introduction/curation`. Catalogue maintenance is described in
{doc}`/developers/documentation`.
