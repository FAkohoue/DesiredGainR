# Security and responsible disclosure

## Supported versions

Security and data-integrity fixes are applied to the current development
version and to the latest released minor version where practical.

## Reporting

Do not open a public issue for a vulnerability that could expose data, execute
untrusted content, or corrupt scientific results. Report it privately to
`akohoue.f@gmail.com` with:

- the affected DesiredGainR and R versions;
- a minimal reproduction using synthetic, non-sensitive data;
- the expected and observed behaviour;
- the likely impact; and
- any proposed mitigation.

Do not include credentials, private repository links, confidential germplasm
identifiers, or real genotype and phenotype records.

## Scientific integrity reports

A defect that silently returns a wrong number is treated with the same urgency
as a conventional vulnerability, because the output of this package feeds
selection decisions. Report privately, by the route above, any defect that:

- returns candidates, coefficients, or gains that are wrong while appearing
  internally consistent;
- misapplies trait direction, units, or row alignment; or
- reports a quantity under a label that does not describe what was computed.

Include the two derived quantities that disagree, where the problem was found
that way, since such defects are rarely visible from a single output.

## External data

DesiredGainR reads only the data supplied to it and invokes no external
executables. Users remain responsible for the provenance of their genotype and
phenotype records and for protecting the directories in which results are
written.
