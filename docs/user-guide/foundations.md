# Foundations

The common conventions are SI units, vacuum wavelengths, and explicit conditions.
For gases, distinguishing a per-particle response from a density-dependent index
is essential. Begin with the {doc}`gas-index example </quickstart>` if you want to
see those conventions in a calculation first.

Read the conventions before comparing optmat with a published value. Then read
about model identity to understand status, qualified intervals, and what to save
with a result. The next chapter applies these ideas to {doc}`linear`.

```{toctree}
:maxdepth: 1

../conventions
../identity
```

Reference interfaces: {mod}`optmat.conditions`, {mod}`optmat.identity`, and
{mod}`optmat.exceptions`. For the source's own units and conditions, see
{doc}`/introduction/curation` and {doc}`/data`.
