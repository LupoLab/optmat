# From a paper to a model

The package distinguishes a source, a dataset, and a model. That distinction lets
you inspect a result without confusing what a paper reported with how optmat uses
it. The {doc}`/quickstart` shows how to retrieve a dataset and its citation in Python.

## Sources and datasets

A **source record** describes one paper, compilation, calculation, or database. It
holds the citation, DOI where available, technique, and notes about the source's
conditions and conventions. A **dataset** records one quantity for one material
from that source, with its units, reference state where relevant, uncertainty, and
location in the source.

Source-specific extraction scripts preserve the strings printed in tables. Values
published only in the text are recorded with their page. Corrections are separate
errata: each keeps the printed value, the correction, the reason, and the evidence.
The loader applies the recorded units and corrections when the data are read.

The records distinguish original data, reproduced values, checks, and exclusions.
This matters when several papers quote the same measurement: they are not
independent evidence. Some records contain formulas or notes without extracted
datasets. Likewise, the presence of data for a material does not mean that optmat
provides a recommended model for every property of that material.

Browse {doc}`/the-data/sources` for the recorded citations and contributions.
{doc}`/data` describes the record format and conversions;
{doc}`/the-data/quantities` defines the vocabulary. {mod}`optmat.datasets` is the Python interface to the records.

## Models and their support

A model may implement a published formula, interpolate a curated table, or fit
selected datasets. A fit records its data selection and uncertainty treatment;
measurements held out of the fit can test its predictions. The model's description
explains the physical assumptions and which data support it.

Use the {doc}`functionality overview <overview>` to find the appropriate model
family. The {doc}`/user-guide/index` explains its conventions, supported ranges,
and validation. Model {doc}`identity, status, and validity </identity>` make these
choices inspectable from a calculation.

Statistical uncertainty, sensitivity to individual datasets, and agreement with
independent measurements answer different questions. Read the relevant accuracy
chapter before treating a numerical estimate as the uncertainty of your application:
{doc}`/accuracy` for gas dispersion or {doc}`/nonlinear_accuracy` for the third-order
response, for example.

## The review register

The {doc}`review register </register>` explains what was checked, corrections made, decisions about
fits, and questions still open. References such as `TC-92` point to individual
entries. The source record describes the data; the register preserves the evidence
and reasoning behind their treatment.

You can examine the {doc}`gallery </auto_gallery/index>` to see data, fits, and
validation together. To contribute a source or reproduce a fit, follow the
{doc}`/developers/index`; to use the existing models, continue with the
{doc}`/user-guide/index`.
