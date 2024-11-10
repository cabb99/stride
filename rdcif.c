#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stride.h"


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
        // Create a new chain
        if (*Cn >= MAX_CHAIN)
        {
            fprintf(stderr, "Exceeded maximum number of chains (%d)\n", MAX_CHAIN);
            return FAILURE;
        }
        c = (CHAIN *)malloc(sizeof(CHAIN));
        if (!c)
        {
            fprintf(stderr, "Memory allocation failed for CHAIN\n");
            return FAILURE;
        }
        memset(c, 0, sizeof(CHAIN));
        c->Id = chain_id;
        c->NRes = 0;
        c->Rsd = (RESIDUE **)malloc(MAX_RES * sizeof(RESIDUE *));
        if (!c->Rsd)
        {
            fprintf(stderr, "Memory allocation failed for CHAIN->Rsd\n");
            free(c);
            return FAILURE;
        }
        memset(c->Rsd, 0, MAX_RES * sizeof(RESIDUE *));
        c->Valid = YES;
        Chain[*Cn] = c;
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
        return FAILURE;
    }

    return SUCCESS;
}



int ReadCIFFile(CHAIN **Chain, int *Cn, COMMAND *Cmd)
{
    printf("%s:%d\n", __FILE__, __LINE__);
    int ChainCnt, InfoCnt, i;
    enum METHOD Method = XRay;
    BOOLEAN First_ATOM, Published = YES, DsspAssigned = NO;
    float Resolution = 0.0;
    FILE *cif;
    BUFFER Buffer;
    char *Info[MAX_INFO], PdbIdent[5];
    RESIDUE *r;
    CHAIN *c;

    *Cn = 0;
    InfoCnt = 0;
    strcpy(PdbIdent, "~~~~");

    if (!(cif = fopen(Cmd->InputFile, "r")))
        return (FAILURE);

    First_ATOM = YES;

    // Variables for parsing
    int in_atom_site_loop = 0;
    char *atom_site_fields[100];
    int num_fields = 0;
    int field_indices[100];
    int num_atoms = 0;

    // Indices for specific atom_site fields
    int idx_group_PDB = -1;
    int idx_id = -1;
    int idx_type_symbol = -1;
    int idx_label_atom_id = -1;
    int idx_label_comp_id = -1;
    int idx_label_asym_id = -1;
    int idx_label_entity_id = -1;
    int idx_label_seq_id = -1;
    int idx_Cartn_x = -1;
    int idx_Cartn_y = -1;
    int idx_Cartn_z = -1;
    int idx_occupancy = -1;
    int idx_B_iso_or_equiv = -1;
    int idx_label_alt_id = -1;
    printf("%s:%d\n", __FILE__, __LINE__);
    // Read the CIF file line by line
    while (fgets(Buffer, BUFSZ, cif))
    {
        if (strncmp(Buffer, "data_", 5) == 0)
        {
            // Process data block
            strncpy(PdbIdent, Buffer + 5, 4);
            PdbIdent[4] = '\0';
        }
        else if (strncmp(Buffer, "#", 1) == 0)
        {
            // Comment or separator, skip
            continue;
        }
        else if (strncmp(Buffer, "loop_", 5) == 0)
        {
            // Start of a loop
            in_atom_site_loop = 0;
            num_fields = 0;

            // Read the data item names
            while (fgets(Buffer, BUFSZ, cif))
            {
                if (Buffer[0] == '_')
                {
                    if (strncmp(Buffer, "_atom_site.", 11) == 0)
                    {
                        atom_site_fields[num_fields] = strdup(Buffer);
                        atom_site_fields[num_fields][strcspn(atom_site_fields[num_fields], " \t\n")] = '\0'; // Remove trailing whitespace
                        num_fields++;
                        in_atom_site_loop = 1;
                    }
                    else
                    {
                        // Not an _atom_site field
                        if (in_atom_site_loop)
                            break; // End of _atom_site fields
                    }
                }
                else
                {
                    // Not a data item name
                    break;
                }
            }
            printf("%s:%d\n", __FILE__, __LINE__);
            if (in_atom_site_loop)
            {
                // Identify indices of required fields
                for (i = 0; i < num_fields; i++)
                {
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

                // Now Buffer contains the first data line
                do
                {
                    // Process Buffer
                    // Split Buffer into tokens
                    char *tokens[100];
                    int num_tokens = 0;
                    char *p = strtok(Buffer, " \t\n");
                    while (p != NULL)
                    {
                        tokens[num_tokens++] = p;
                        p = strtok(NULL, " \t\n");
                    }

                    if (num_tokens == num_fields)
                    {
                        // Map tokens to fields
                        // Create an ATOM record and add it to the Chain
                        if (!Process_CIF_ATOM(tokens, num_fields, idx_group_PDB, idx_id, idx_type_symbol,
                                              idx_label_atom_id, idx_label_comp_id, idx_label_asym_id, idx_label_entity_id,
                                              idx_label_seq_id, idx_Cartn_x, idx_Cartn_y, idx_Cartn_z, idx_occupancy,
                                              idx_B_iso_or_equiv, idx_label_alt_id, Chain, Cn, &First_ATOM, Cmd))
                        {
                            return (FAILURE);
                        }
                        num_atoms++;
                    }
                    else
                    {
                        // The number of tokens does not match the number of fields
                        // Probably end of data
                        break;
                    }
                } while (fgets(Buffer, BUFSZ, cif) && Buffer[0] != '#' && Buffer[0] != '_');
            }
        }
        // Process other data items as needed
        else if (strncmp(Buffer, "_struct.title", 13) == 0)
        {
            // Process title
            Info[InfoCnt] = (char *)ckalloc(BUFSZ * sizeof(char));
            strcpy(Info[InfoCnt], "HDR  ");
            char *title = Buffer + 13; // Skip '_struct.title'
            strcat(Info[InfoCnt++], title);
        }
        else if (strncmp(Buffer, "_citation.author", 16) == 0)
        {
            // Process author
            Info[InfoCnt] = (char *)ckalloc(BUFSZ * sizeof(char));
            strcpy(Info[InfoCnt], "AUT  ");
            char *author = Buffer + 16; // Skip '_citation.author'
            strcat(Info[InfoCnt++], author);
        }
        // Continue processing other data items
    }
    printf("%s:%d\n", __FILE__, __LINE__);
    fclose(cif);
    printf("%s:%d\n", __FILE__, __LINE__);
    // Process the rest of the data and fill the Chain structures
    for (ChainCnt = 0; ChainCnt < *Cn; ChainCnt++)
    {
        printf("%s:%d\n", __FILE__, __LINE__);
        c = Chain[ChainCnt];
        printf("%s:%d\n", __FILE__, __LINE__);
        printf("Debug: Calling FindAtom with chain %c, NRes %d, atom name %s\n", c->Id, c->NRes, "CA");
        if (c->NRes != 0 && !FindAtom(c, c->NRes-1, "CA", &i))
            c->NRes--;
        printf("%s:%d\n", __FILE__, __LINE__);
        printf("Debug: Chain File: %s\n", c->File);
        printf("Debug: Command Input File: %s\n", Cmd->InputFile);
        if (c->File && Cmd->InputFile)
        {
            printf("Debug: Both Chain File and Command Input File exist.\n");
        }
        else
        {
            printf("Debug: One or both of Chain File and Command Input File do not exist.\n");
        }

        if (c->File == NULL) {
            c->File = (char *)malloc(BUFSZ * sizeof(char));
            if (!c->File) {
                fprintf(stderr, "Memory allocation failed for c->File\n");
                return FAILURE;
            }
        }
        strcpy(c->File, Cmd->InputFile);
        printf("%s:%d\n", __FILE__, __LINE__);
        strcpy(c->PdbIdent, PdbIdent);
        // if (c->NRes != 0)
        //     c->NRes++;
        if (c->NSheet != -1)
            c->NSheet++;
        c->Resolution = Resolution;
        c->Method = Method;
        c->Published = Published;
        c->DsspAssigned = DsspAssigned;
        c->NInfo = InfoCnt;
        printf("%s:%d\n", __FILE__, __LINE__);
        if (c->Info == NULL) {
            c->Info = (char **)malloc(MAX_INFO * sizeof(char *));
            if (!c->Info) {
                fprintf(stderr, "Memory allocation failed for c->Info\n");
                return FAILURE;
            }
            memset(c->Info, 0, MAX_INFO * sizeof(char *)); // Optional, for safety
        }
        
        for (i = 0; i < InfoCnt; i++)
        {
            printf("%s:%d\n", __FILE__, __LINE__);
            c->Info[i] = (char *)ckalloc(BUFSZ * sizeof(char));
            printf("%s:%d\n", __FILE__, __LINE__);
            strcpy(c->Info[i], Info[i]);
            printf("%s:%d\n", __FILE__, __LINE__);
            c->Info[i][71] = '\0';
        }
        printf("%s:%d\n", __FILE__, __LINE__);
            for (i = 0; i < InfoCnt; i++)
        free(Info[i]);
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

        for (i = 0; i < c->NRes; i++)
        {
            // printf("Debug: Processing residue %d in chain %c\n", i, c->Id);
            r = c->Rsd[i];
            if (r == NULL) {
            fprintf(stderr, "Error: Residue %d in chain %c is NULL\n", i, c->Id);
            return FAILURE;
            }
            r->Inv = (INVOLVED *)ckalloc(sizeof(INVOLVED));
            if (r->Inv == NULL) {
            fprintf(stderr, "Error: Memory allocation failed for r->Inv\n");
            return FAILURE;
            }
            r->Prop = (PROPERTY *)ckalloc(sizeof(PROPERTY));
            if (r->Prop == NULL) {
            fprintf(stderr, "Error: Memory allocation failed for r->Prop\n");
            return FAILURE;
            }
            r->Inv->NBondDnr = 0;
            r->Inv->NBondAcc = 0;
            r->Inv->InterchainHBonds = NO;
            r->Prop->Asn = 'C';
            r->Prop->PdbAsn = 'C';
            r->Prop->DsspAsn = 'C';
            r->Prop->Solv = 0.0;
            r->Prop->Phi = 360.0;
            r->Prop->Psi = 360.0;
            // printf("Debug: Finished processing residue %d in chain %c\n", i, c->Id);
        }
    }
    printf("%s:%d\n", __FILE__, __LINE__);
    // for (i = 0; i < InfoCnt; i++)
    //     free(Info[i]);
    //         printf("Chain ID: %c\n", c->Id);
    //     printf("Number of Residues: %d\n", c->NRes);
    //     printf("Number of Atoms: %d\n", c->NAtom);
    //     printf("Resolution: %.2f\n", c->Resolution);
    //     printf("Method: %d\n", c->Method);
    //     printf("Published: %d\n", c->Published);
    //     printf("DSSP Assigned: %d\n", c->DsspAssigned);
    //     printf("Number of Info Records: %d\n", c->NInfo);
    //     for (i = 0; i < c->NInfo; i++) {
    //         printf("Info[%d]: %s\n", i, c->Info[i]);
    //     }
    //     printf("Residue Number | Residue Name | Atom Name | X Coordinate | Y Coordinate | Z Coordinate | Occupancy | B Factor\n");
    //     printf("---------------------------------------------------------------------------------------------\n");
    //     for (i = 0; i < c->NRes; i++) {
    //         r = c->Rsd[i];
    //         for (int j = 0; j < r->NAtom; j++) {
    //             printf("%-14s | %-12s | %-9s | %-12.3f | %-12.3f | %-12.3f | %-9.2f | %-8.2f\n",
    //                    r->PDB_ResNumb, r->ResType, r->AtomType[j], r->Coord[j][0], r->Coord[j][1], r->Coord[j][2], r->Occupancy[j], r->TempFactor[j]);
    //         }
    //     }
    // printf("%s:%d\n", __FILE__, __LINE__);
    return (SUCCESS);
}