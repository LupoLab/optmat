# Enhancement proposals

An optmat Enhancement Proposal (OEP) is the design document for one major change. That
might be:
- a new medium;
- a new family of properties or models;
- a new kind of curated value;
- a changed default or convention;
- a redesigned interface.

An OEP sets out the motivation, the theory with its references, the primary sources it
needs, the implementation options and how the result will be verified. Once accepted,
it is implemented in phases, and each phase appends its report to the OEP. Finished OEPs
are archived and linked from the manual pages they produced.

[OEP 0](oep-0000-process.md) describes the process and says when an OEP is needed. New
OEPs start from [the template](template.md). Every OEP is reviewed by at least one human
and one frontier AI model before it is accepted.

The {doc}`roadmap <../roadmap>` lists, in order:
1. accepted and in-progress OEPs;
2. drafts under discussion;
3. informal work not fully covered by a proposal.

## Open

Proposals that are active, in draft, accepted, in progress or deferred.

| OEP | Title | Status | Type |
|---|---|---|---|
| [0](oep-0000-process.md) | Purpose and process | Active | Process |
| [7](oep-0007-strong-field-dissociation.md) | Strong-field dissociation of small molecules: the products a pulse leaves | Draft | Feature |
| [8](oep-0008-bound-two-photon-response.md) | The bound two-photon resonances of the noble gases for propagation | Draft | Feature |

## Archived

Proposals that are final, rejected, withdrawn or superseded: see the
[archive](archive/index.md).

| OEP | Title | Status | Type |
|---|---|---|---|
| [1](archive/oep-0001-propagation-interface.md) | The propagation-code interface: identity, overrides and the data for solprop's SEP 3 | Final | Feature |
| [2](archive/oep-0002-compact-infrared.md) | Compact infrared terms for the molecular gases | Final | Feature |
| [3](archive/oep-0003-infrared-windows.md) | The molecular gases between their infrared bands | Final | Feature |
| [4](archive/oep-0004-mixture-index.md) | The refractive index of gas mixtures at finite density | Final | Feature |
| [5](archive/oep-0005-hyperpolarizability-support.md) | A public record of where the hyperpolarizability fits are supported | Final | Feature |
| [6](archive/oep-0006-two-photon-ionisation.md) | Two-photon ionisation and absorption of the noble gases | Final | Feature |

```{toctree}
:hidden:
:maxdepth: 1

oep-0000-process
oep-0007-strong-field-dissociation
oep-0008-bound-two-photon-response
template
archive/index
```
