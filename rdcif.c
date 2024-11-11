#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "stride.h"


#define MAX_CIF_TOKENS 100

int cif_tokenize(const char *line, char **tokens, int max_tokens) {
    int num_tokens = 0;
    const char *p = line;
    while (*p && num_tokens < max_tokens) {
        // Skip leading whitespace
        while (*p && isspace(*p))
            p++;

        if (*p == '\0' || *p == '#')
            break;

        char *token;
        if (*p == '\'' || *p == '"') {
            char quote = *p++;
            const char *start = p;
            while (*p && *p != quote)
                p++;
            size_t len = p - start;
            token = (char *)ckalloc(len + 1);
            strncpy(token, start, len);
            token[len] = '\0';
            if (*p == quote)
                p++;
        } else if (*p == ';') {
            // Handle semicolon-delimited strings (multi-line)
            p++; // Skip the semicolon
            const char *start = p;
            while (*p && !(*p == ';' && (*(p - 1) == '\n' || *(p - 1) == '\r')))
                p++;
            size_t len = p - start;
            token = (char *)ckalloc(len + 1);
            strncpy(token, start, len);
            token[len] = '\0';
            if (*p == ';')
                p++;
        } else {
            const char *start = p;
            while (*p && !isspace(*p))
                p++;
            size_t len = p - start;
            token = (char *)ckalloc(len + 1);
            strncpy(token, start, len);
            token[len] = '\0';
        }
        tokens[num_tokens++] = token;
    }
    return num_tokens;
}

int IsAminoAcid(char *residue_name)
{
    const char *amino_acids[] = {
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
        "ASX", "GLX", "SEC", "PYL" // Include any others as needed
    };
    int num_amino_acids = sizeof(amino_acids) / sizeof(amino_acids[0]);
    int i;
    for (i = 0; i < num_amino_acids; i++)
    {
        if (strcmp(residue_name, amino_acids[i]) == 0)
            return 1;
    }
    return 0;
}

int AddAtomToChain(CHAIN **Chain, int *Cn, char *atom_name, char alt_loc,
                   char *residue_name, char *residue_number, char chain_id,
                   float x, float y, float z, float occupancy, float b_factor,
                   char *element, BOOLEAN *First_ATOM)
{
    CHAIN *c = NULL;
    RESIDUE *r = NULL;
    int i;

    // Find or create the chain
    for (i = 0; i < *Cn; i++)
    {
        if (Chain[i]->Id == chain_id)
        {
            c = Chain[i];
            break;
        }
    }

    if (c == NULL)
    {
        // Create a new chain using InitChain
        if (*Cn >= MAX_CHAIN)
        {
            fprintf(stderr, "Exceeded maximum number of chains (%d)\n", MAX_CHAIN);
            return FAILURE;
        }
        InitChain(&Chain[*Cn]);
        c = Chain[*Cn];
        c->Id = chain_id;
        (*Cn)++;
    }

    // Find or create the residue
    for (i = 0; i < c->NRes; i++)
    {
        if (strcmp(c->Rsd[i]->PDB_ResNumb, residue_number) == 0)
        {
            r = c->Rsd[i];
            break;
        }
    }

    if (r == NULL)
    {
        // Create a new residue
        if (c->NRes >= MAX_RES)
        {
            fprintf(stderr, "Exceeded maximum number of residues (%d) in chain %c\n", MAX_RES, chain_id);
            return FAILURE;
        }
        r = (RESIDUE *)malloc(sizeof(RESIDUE));
        if (!r)
        {
            fprintf(stderr, "Memory allocation failed for RESIDUE\n");
            return FAILURE;
        }
        memset(r, 0, sizeof(RESIDUE));
        strncpy(r->ResType, residue_name, RES_FIELD - 1);
        r->ResType[RES_FIELD - 1] = '\0';
        strncpy(r->PDB_ResNumb, residue_number, RES_FIELD - 1);
        r->PDB_ResNumb[RES_FIELD - 1] = '\0';
        r->NAtom = 0;
        c->Rsd[c->NRes] = r;
        c->NRes++;
    }

    // Add the atom to the residue
    if (r->NAtom >= MAX_AT_IN_RES)
    {
        fprintf(stderr, "Exceeded maximum number of atoms (%d) in residue %s\n", MAX_AT_IN_RES, residue_number);
        return FAILURE;
    }

    strncpy(r->AtomType[r->NAtom], atom_name, AT_FIELD - 1);
    r->AtomType[r->NAtom][AT_FIELD - 1] = '\0';

    r->Coord[r->NAtom][0] = x;
    r->Coord[r->NAtom][1] = y;
    r->Coord[r->NAtom][2] = z;
    r->Occupancy[r->NAtom] = occupancy;
    r->TempFactor[r->NAtom] = b_factor;

    // Note: The element symbol and alt_loc are not stored in the given structures.
    // If necessary, you might need to adjust your structures to include them.

    r->NAtom++;

    // Update chain's NAtom
    c->NAtom++;

    return SUCCESS;
}


int Process_CIF_ATOM(char **tokens, int num_fields,
                     int idx_group_PDB, int idx_id, int idx_type_symbol,
                     int idx_label_atom_id, int idx_label_comp_id,
                     int idx_label_asym_id, int idx_label_entity_id,
                     int idx_label_seq_id, int idx_Cartn_x,
                     int idx_Cartn_y, int idx_Cartn_z, int idx_occupancy,
                     int idx_B_iso_or_equiv, int idx_label_alt_id,
                     CHAIN **Chain, int *Cn, BOOLEAN *First_ATOM, COMMAND *Cmd)
{
    float x, y, z, occupancy, b_factor;
    char element[3], atom_name[5], residue_name[4], chain_id;
    char residue_number[RES_FIELD]; // PDB residue number as string
    char alt_loc;

    // Parse required fields
    // Coordinates are mandatory
    if (idx_Cartn_x != -1 && idx_Cartn_y != -1 && idx_Cartn_z != -1)
    {
        x = atof(tokens[idx_Cartn_x]);
        y = atof(tokens[idx_Cartn_y]);
        z = atof(tokens[idx_Cartn_z]);
    }
    else
    {
        // Coordinates are mandatory
        printf("Error: Missing coordinates in ATOM record\n");
        printf("Debug: %d %d %d\n", idx_Cartn_x, idx_Cartn_y, idx_Cartn_z);
        return FAILURE;
    }

    occupancy = (idx_occupancy != -1) ? atof(tokens[idx_occupancy]) : 1.0;
    b_factor = (idx_B_iso_or_equiv != -1) ? atof(tokens[idx_B_iso_or_equiv]) : 0.0;

    if (idx_type_symbol != -1)
    {
        strncpy(element, tokens[idx_type_symbol], 2);
        element[2] = '\0';
    }
    else
        strcpy(element, "  ");

    if (idx_label_atom_id != -1)
    {
        strncpy(atom_name, tokens[idx_label_atom_id], 4);
        atom_name[4] = '\0';
    }
    else
        strcpy(atom_name, "    ");

    if (idx_label_comp_id != -1)
    {
        strncpy(residue_name, tokens[idx_label_comp_id], 3);
        residue_name[3] = '\0';
    }
    else
        strcpy(residue_name, "UNK");

    if (!IsAminoAcid(residue_name))
    {
        // Skip non-amino acid residues
        return SUCCESS;
    }

    if (idx_label_seq_id != -1)
    {
        printf("Debug: Processing residue number: %s\n", tokens[idx_label_seq_id]);
        strncpy(residue_number, tokens[idx_label_seq_id], RES_FIELD - 1);
        residue_number[RES_FIELD - 1] = '\0';
    }
    else
        strcpy(residue_number, "   "); // Default or error

    chain_id = (idx_label_asym_id != -1) ? tokens[idx_label_asym_id][0] : ' ';
    alt_loc = (idx_label_alt_id != -1) ? tokens[idx_label_alt_id][0] : ' ';

    // Now, process the atom and add it to the data structures
    if (!AddAtomToChain(Chain, Cn, atom_name, alt_loc, residue_name, residue_number,
                        chain_id, x, y, z, occupancy, b_factor, element, First_ATOM))
    {
        printf("Error: Failed to add atom to chain\n");
        printf("Debug: %s %c %s %s %c %f %f %f %f %f %s\n",
               atom_name, alt_loc, residue_name, residue_number, chain_id,
               x, y, z, occupancy, b_factor, element);
        return FAILURE;
    }

    return SUCCESS;
}



int ReadCIFFile(CHAIN **Chain, int *Cn, COMMAND *Cmd)
{
    int ChainCnt, InfoCnt = 0, i;
    enum METHOD Method = XRay;
    BOOLEAN First_ATOM = YES, Published = YES, DsspAssigned = NO;
    float Resolution = 0.0;
    FILE *cif;
    BUFFER Buffer;
    char *Info[MAX_INFO], PdbIdent[5] = "~~~~";
    RESIDUE *r;
    CHAIN *c;

    *Cn = 0;
    
    if (!(cif = fopen(Cmd->InputFile, "r"))) {
        fprintf(stderr, "Failed to open CIF file: %s\n", Cmd->InputFile);
        return FAILURE;
    }

    // Variables for parsing
    int in_atom_site_loop = 0;
    char *atom_site_fields[MAX_CIF_TOKENS];
    int num_fields = 0;

    // Indices for specific atom_site fields
    int idx_group_PDB = -1, idx_id = -1, idx_type_symbol = -1, idx_label_atom_id = -1;
    int idx_label_comp_id = -1, idx_label_asym_id = -1, idx_label_entity_id = -1;
    int idx_label_seq_id = -1, idx_Cartn_x = -1, idx_Cartn_y = -1, idx_Cartn_z = -1;
    int idx_occupancy = -1, idx_B_iso_or_equiv = -1, idx_label_alt_id = -1;

    
    while (fgets(Buffer, BUFSZ, cif)) {
        // Remove trailing newline
        Buffer[strcspn(Buffer, "\n")] = '\0';

        if (strncmp(Buffer, "data_", 5) == 0) {
            strncpy(PdbIdent, Buffer + 5, 4);
            PdbIdent[4] = '\0';
        } else if (strncmp(Buffer, "_struct.title", 13) == 0) {
            Info[InfoCnt] = (char *)ckalloc(BUFSZ * sizeof(char));
            strcpy(Info[InfoCnt], "HDR  ");
            strcat(Info[InfoCnt++], Buffer + 14);
        } else if (strncmp(Buffer, "_citation.author", 16) == 0) {
            Info[InfoCnt] = (char *)ckalloc(BUFSZ * sizeof(char));
            strcpy(Info[InfoCnt], "AUT  ");
            strcat(Info[InfoCnt++], Buffer + 17);
        } else if (strncmp(Buffer, "loop_", 5) == 0) {
            // Start of a loop
            num_fields = 0;
            in_atom_site_loop = 0;

            // Read data item names
            while (fgets(Buffer, BUFSZ, cif)) {
                Buffer[strcspn(Buffer, "\n")] = '\0';

                // Skip empty lines
                if (strlen(Buffer) == 0)
                    continue;

                // Check if line starts with '_', indicating data item names
                if (Buffer[0] == '_') {
                    // Tokenize the line in case multiple data item names are on the same line
                    char *line = Buffer;
                    char *token;
                    while ((token = strtok_r(line, " \t", &line))) {
                        if (strncmp(token, "_atom_site.", 11) == 0) {
                            atom_site_fields[num_fields] = strdup(token);
                            num_fields++;
                            in_atom_site_loop = 1;
                        } else if (in_atom_site_loop) {
                            // Finished reading _atom_site data names
                            break;
                        }
                    }
                } else if (in_atom_site_loop) {
                    // We've finished reading data item names and are now at data lines
                    break;
                }
            }

            if (in_atom_site_loop && num_fields > 0) {
                // Map field names to indices
                for (i = 0; i < num_fields; i++) {
                    if (strcmp(atom_site_fields[i], "_atom_site.group_PDB") == 0)
                        idx_group_PDB = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.id") == 0)
                        idx_id = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.type_symbol") == 0)
                        idx_type_symbol = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.label_atom_id") == 0)
                        idx_label_atom_id = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.label_comp_id") == 0)
                        idx_label_comp_id = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.label_asym_id") == 0)
                        idx_label_asym_id = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.label_entity_id") == 0)
                        idx_label_entity_id = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.label_seq_id") == 0)
                        idx_label_seq_id = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.Cartn_x") == 0)
                        idx_Cartn_x = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.Cartn_y") == 0)
                        idx_Cartn_y = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.Cartn_z") == 0)
                        idx_Cartn_z = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.occupancy") == 0)
                        idx_occupancy = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.B_iso_or_equiv") == 0)
                        idx_B_iso_or_equiv = i;
                    else if (strcmp(atom_site_fields[i], "_atom_site.label_alt_id") == 0)
                        idx_label_alt_id = i;
                }

                // Process data lines
                do {
                    Buffer[strcspn(Buffer, "\n")] = '\0';
                    // Skip empty lines and comments
                    if (Buffer[0] == '#' || strlen(Buffer) == 0)
                        continue;
                    if (Buffer[0] == '_') {
                        // Start of new data item names or loop, break
                        break;
                    }

                    char *tokens[MAX_CIF_TOKENS];
                    int num_tokens = cif_tokenize(Buffer, tokens, MAX_CIF_TOKENS);
                    if (num_tokens == num_fields) {
                        if (!Process_CIF_ATOM(tokens, num_fields, idx_group_PDB, idx_id, idx_type_symbol,
                                              idx_label_atom_id, idx_label_comp_id, idx_label_asym_id, idx_label_entity_id,
                                              idx_label_seq_id, idx_Cartn_x, idx_Cartn_y, idx_Cartn_z, idx_occupancy,
                                              idx_B_iso_or_equiv, idx_label_alt_id, Chain, Cn, &First_ATOM, Cmd)) {
                            fclose(cif);
                            fprintf(stderr, "Error: Failed to process ATOM record\n");
                            return FAILURE;
                        }
                    } else {
                        // Mismatch in number of tokens and fields
                        fprintf(stderr, "Mismatch in tokens and fields at line: %s\n", Buffer);
                    }

                    // Free tokens
                    for (i = 0; i < num_tokens; i++)
                        free(tokens[i]);

                } while (fgets(Buffer, BUFSZ, cif));

                // Free field names
                for (i = 0; i < num_fields; i++)
                    free(atom_site_fields[i]);
            }
        }
    }

    fclose(cif);
    printf("%s:%d\n", __FILE__, __LINE__);

    // Process the rest of the data and fill the Chain structures
    for (ChainCnt = 0; ChainCnt < *Cn; ChainCnt++) {
        printf("%s:%d\n", __FILE__, __LINE__);
        c = Chain[ChainCnt];
        printf("%s:%d\n", __FILE__, __LINE__);
        if (c->NRes != 0 && !FindAtom(c, c->NRes - 1, "CA", &i))
            c->NRes--;
        printf("%s:%d\n", __FILE__, __LINE__);
        printf("Cmd->InputFile: %s\n", Cmd->InputFile);
        c->File = (char *)ckalloc(BUFSZ * sizeof(char));
        if (!c->File) {
            fprintf(stderr, "Memory allocation failed for c->File\n");
            return FAILURE;
        }
        strcpy(c->File, Cmd->InputFile);
        printf("%s:%d\n", __FILE__, __LINE__);
        strcpy(c->PdbIdent, PdbIdent);
        printf("%s:%d\n", __FILE__, __LINE__);
        printf("%s:%d:%d\n", __FILE__, __LINE__, c->NRes);
        // if (c->NRes != 0)
        //     c->NRes++;
        printf("%s:%d:%d\n", __FILE__, __LINE__, c->NRes);
        printf("%s:%d\n", __FILE__, __LINE__);
        if (c->NSheet != -1)
            c->NSheet++;
        printf("%s:%d\n", __FILE__, __LINE__);
        c->Resolution = Resolution;
        c->Method = Method;
        c->Published = Published;
        c->DsspAssigned = DsspAssigned;
        c->NInfo = InfoCnt;
        printf("%s:%d\n", __FILE__, __LINE__);
        for (i = 0; i < InfoCnt; i++) {
            c->Info[i] = (char *)ckalloc(BUFSZ * sizeof(char));
            strcpy(c->Info[i], Info[i]);
            c->Info[i][71] = '\0';
        }
        printf("%s:%d\n", __FILE__, __LINE__);
        for (i = 0; i < c->NRes; i++) {
            printf("%s:%d\n", __FILE__, __LINE__);
            r = c->Rsd[i];
            printf("%s:%d:%d/%d\n", __FILE__, __LINE__, i, c->NRes);
            r->Inv = (INVOLVED *)ckalloc(sizeof(INVOLVED));
            printf("%s:%d\n", __FILE__, __LINE__);
            r->Prop = (PROPERTY *)ckalloc(sizeof(PROPERTY));
            printf("%s:%d\n", __FILE__, __LINE__);
            r->Inv->NBondDnr = 0;
            r->Inv->NBondAcc = 0;
            r->Inv->InterchainHBonds = NO;
            r->Prop->Asn = 'C';
            r->Prop->PdbAsn = 'C';
            r->Prop->DsspAsn = 'C';
            r->Prop->Solv = 0.0;
            r->Prop->Phi = 360.0;
            r->Prop->Psi = 360.0;
        }
        printf("%s:%d\n", __FILE__, __LINE__);
        printf("Chain ID: %c\n", c->Id);
        printf("Number of Residues: %d\n", c->NRes);
        printf("Number of Atoms: %d\n", c->NAtom);
        printf("Resolution: %.2f\n", c->Resolution);
        printf("Method: %d\n", c->Method);
        printf("Published: %d\n", c->Published);
        printf("DSSP Assigned: %d\n", c->DsspAssigned);
        printf("Number of Info Records: %d\n", c->NInfo);
        for (i = 0; i < c->NInfo; i++) {
            printf("Info[%d]: %s\n", i, c->Info[i]);
        }
        printf("Residue Number | Residue Name | Atom Name | X Coordinate | Y Coordinate | Z Coordinate | Occupancy | B Factor\n");
        printf("---------------------------------------------------------------------------------------------\n");
        for (i = 0; i < c->NRes; i++) {
            r = c->Rsd[i];
            for (int j = 0; j < r->NAtom; j++) {
                printf("%-14s | %-12s | %-9s | %-12.3f | %-12.3f | %-12.3f | %-9.2f | %-8.2f\n",
                        r->PDB_ResNumb, r->ResType, r->AtomType[j], r->Coord[j][0], r->Coord[j][1], r->Coord[j][2], r->Occupancy[j], r->TempFactor[j]);
            }
        }
        printf("%s:%d\n", __FILE__, __LINE__);
    }
    for (i = 0; i < InfoCnt; i++)
        free(Info[i]);
    printf("%s:%d\n", __FILE__, __LINE__);
    return SUCCESS;
}