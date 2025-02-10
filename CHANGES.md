Modified STRIDE – with CIF File Support
===============================================

This version of STRIDE incorporates modifications to add support for 
reading and processing CIF files in addition to the standard PDB format. Below 
is a summary of the main changes:

   - Added a new source file, **rdcif.c**, which contains routines to parse CIF files.
   - Updated the main routine (in **stride.c**) to detect the file extension and call 
     the appropriate file-reading function (i.e. `ReadPDBFile` for “.pdb” files and 
     `ReadCIFFile` for “.cif” files).
   - Modified the **Makefile** to include the new source file `rdcif.c` in the build.
   - Changed functions that previously assumed chain identifiers as single 
     characters. Now, chain identifiers are treated as strings (e.g. changes in 
     `chkchain.c`, `contact_map.c`, `contact_order.c`, `dssp.c`, etc.).
   - Added utility functions such as `ChainInList` in **strutil.c** and a new version 
     of `SpaceToDash` that returns a string.
   - Updated command-line option parsing (in **stride.c**) to allow chain selection using 
     both the old (each letter is one chain) and new (comma-separated, multi–letter IDs) 
     syntaxes.
   - Added helper functions `parseChainList` and `parseOldChainList` for backward compatibility
   - Revised printing formats (changed format specifiers from `%c` to `%s` for chain IDs) 
   - Other small corrections (e.g. initialization changes in `initchn.c` and minor 
     adjustments in string handling)