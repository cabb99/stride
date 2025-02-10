# Modified STRIDE with CIF Support

This repository contains a modified version of the STRIDE program, originally developed by Dmitrij Frishman and Patrick Argos at the European Molecular Biology Laboratory (EMBL) in Heidelberg, Germany, for secondary structure prediction of protein structures. In this version, the code has been modified to support reading and processing mmCIF files in addition to the standard PDB format.

## Installation

This version of STRIDE is distributed as source code with a Makefile for compilation.

### Compilation

Clone the repository and run the following command in the project directory:

```bash
make
```

This will compile the source files and produce an executable (typically named `stride`).

### Usage

The Method is presented in detail in:

Frishman D, Argos P. Knowledge-Based Protein Secondary Structure Assignment Proteins: Structure, Function, and Genetics 23:566-579 (1995)

For a detailed information on the output refer to https://webclu.bio.wzw.tum.de/stride/ or see the [stride documentation](https://webclu.bio.wzw.tum.de/stride/stride.doc)

For a list of options, run:

```bash
./stride --help
```

---

## Changelog

For a detailed summary of modifications, please see the [CHANGES](CHANGES.md) file included in this repository.

---

## License

This software is provided under the original STRIDE license. For full license details, please see the [LICENSE](LICENSE) file.
