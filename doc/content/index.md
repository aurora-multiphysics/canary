!config navigation breadcrumbs=False scrollspy=False

# Canary

Canary is a collection of tested UKAEA MOOSE demonstration problems, examples, and their
documentation. Every problem here runs with the base [MOOSE](https://mooseframework.inl.gov)
framework alone: Canary ships no application executable of its own, and both the tests and this
documentation are generated using the `moose_test` application found in `MOOSE_DIR`.

## Contents

- [Demonstration problems](problems/index.md)

## Building this documentation id=building

```bash
export MOOSE_DIR=/path/to/moose
cd doc
./moosedocs.py build --serve
```
