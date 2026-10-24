# Adding sources and materials

Read {doc}`/data`, {doc}`/the-data/quantities` and {doc}`/register` before adding
data. The [working rules](https://github.com/LupoLab/optmat/blob/main/AGENTS.md)
define the required evidence; this page gathers the workflow.

## A source or correction

1. Create `src/optmat/data/sources/<source>.toml`. Record the citation and DOI when
   available, technique, measurement conditions, density and wavelength conventions,
   uncertainties and their levels, and the location of every dataset or formula.
2. Write `tools/extract/<source>.py` to extract the printed tables. Preserve the
   printed strings. Values found only in the text go in inline `points` with their
   page. Reproduced values name their original source and are not fitted beside it.
3. Check the table image wherever extraction warnings or structural ambiguities
   affect values. Prefer printed numbers or supplied data to digitisation. The
   figure-extraction rules and approval requirements are in the working rules
   (TC-699); raster values are never fitted.
4. Express corrections as errata with the printed value, replacement, reason and
   evidence. Re-run extraction to change a data file; never hand-edit it.
5. Add {doc}`register` entries for the transcription, checks, corrections and doubts.
6. Rewrite the catalogue index with `uv run python tools/catalogue_index.py`. The
   index, `src/optmat/data/catalogue_index.json`, tells a search by material,
   quantity or formula which records to read; a test fails until it matches the
   records, after any edit to one. An installed copy whose index is out of date warns
   and reads the records themselves.
7. Run the extraction lane and full {doc}`testing` gate. If the data affect a fit,
   also follow {doc}`fitting`.

## Reproducing extraction

Extraction scripts read the curators' paper library, prepared with
[paperextract](https://github.com/jtravs/paperextract). The source record's `library`
key identifies its directory. The package does not need this private library at
runtime.

```sh
PAPER_LIBRARY=<path> uv run poe test-library
```

The tests re-run the scripts into temporary output and compare with the committed
files. If the paper is unavailable in the library, state that limitation; a maintainer
must reproduce the extraction before the change can be verified. A transcription
must still be checked against its primary source, not reconstructed from memory.

## A material or quantity

Register the material in `src/optmat/data/materials.toml` before adding its data.
For a new gas, add the gas entry, sources, fit configuration and stored fit,
default-model choice with its reason, gallery example, accuracy rows and register
entries. Apply the fitting and verification requirements above.

A new quantity needs a vocabulary definition and documentation before data use it.
A new quantity family, medium or model family follows the
{doc}`/oeps/oep-0000-process` process. Shared functionality belongs in the core only
when a second medium needs it; the core must not import a medium's subpackage.
