# About optmat

optmat is a materials library for reliable simulations of ultrafast and nonlinear
optics. Its purpose is to make the route from a published measurement to a model's
prediction explicit, testable and repeatable.

We work back through the primary literature to resolve transcription errors, interpret
measurement conditions correctly, establish where models apply and systematically
evaluate the available data. optmat provides reference implementations of published
fits alongside new fits based on the best available data and analysis. Each keeps its
provenance, conventions, stated range and the evidence for its accuracy.

## Precise and reliable paper extraction

We use [paperextract](https://github.com/jtravs/paperextract) to turn research papers
into a structured library of tables, equations and figures, retaining the original
papers and links to their pages. This supports precise, reliable extraction by making
the source of a table cell and the evidence needed to check it available together.

optmat's source-specific scripts read those extracted tables and preserve their printed
strings. Table images are checked where warnings or structural ambiguities could affect
the numbers, and numerical checks reproduce the paper's own derived values where
available. Values found only in the text are recorded with their page. Corrections keep
the printed value, the corrected value, the reason and the evidence; they never silently
replace the original. The {doc}`curated data <data>` page describes the records and extraction rules.

## Automated verification and reproducible fits

Executable examples, generated-table checks and numerical references connect the
manual to the implementation. The development gate requires 100 % statement and
branch coverage of the package, with the existing exclusions. Coverage describes
which code ran; assertions and independent evidence establish what was checked.
Technical procedures and verification lanes are in {doc}`/developers/testing`.

Fit configurations identify selected datasets, conventions, uncertainty treatment,
priors and ranges. Stored results include covariance and validation information;
the register explains the decisions. See {doc}`fitting` and {doc}`nonlinear` for the
scientific methods and {doc}`/developers/fitting` for reproduction commands.

## Acknowledgements

[refractiveindex.info](https://refractiveindex.info) is an invaluable resource for the
optics community. We are grateful to its maintainers and contributors, and we do not
seek to replace it. Our goals are different and complementary: optmat focuses on
systematic evaluation of the underlying measurements, physical conventions, reference
implementations and new fits with explicit validation and reproducible provenance.

## Contributing

The repository is [LupoLab/optmat on GitHub](https://github.com/LupoLab/optmat).
Corrections, new data, independent comparisons and reports of unexpected results are
welcome. Include the model, conditions, version and a small reproducing example when
reporting a problem; the {doc}`contributing` page explains how to help.

## Use of AI

AI agents assist with literature curation, implementation, tests and documentation,
working under the direction of John C. Travers. He sets the scientific priorities and
verification standards and decides which data, models and conventions are accepted.
AI assistance makes it practical for a small research team to examine a wide body of
literature, derive models, build reference comparisons and maintain detailed
documentation.

The reliability of the result rests on the evidence and checks recorded with it:

- **Primary sources and explicit provenance.** Each curated value names its source,
  location, units, conditions and stated uncertainty. Reproduced values are identified
  so that they cannot silently count as independent measurements.
- **Checks that can be inspected.** The {doc}`review register <register>` records what was transcribed,
  what it was checked against, corrections with their evidence, fitting decisions and
  open questions. It also describes the division of work between AI agents and the
  human maintainer.
- **Independent numerical references.** Tests compare implementations with published
  numerical results, analytic identities and independent calculations. Model validation
  uses measurements held out of the fit where available.
- **Repeatable calculations.** Extraction scripts and fit configurations keep the
  route from a source to a stored result executable. Automated checks detect changes in
  that route and in the results reported in the manual.

Scientific uncertainty and unresolved disagreements remain part of the record. A
successful software check establishes a specific property of an implementation; the
model's stated uncertainty and validation describe how well it represents the material.
