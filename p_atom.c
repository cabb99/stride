#include "stride.h"

int Process_ATOM(BUFFER Buffer, CHAIN **Chain, int *ChainNumber, 
		 BOOLEAN *First_ATOM, COMMAND *Cmd)
{

  char *Field[MAX_FIELD];
  BUFFER Tmp;
  int CC, NR, NA;
  static char LastRes[MAX_CHAIN][RES_FIELD];
  RESIDUE *r;
  
  if( Cmd->NActive && !ChInStr(Cmd->Active,SpaceToDashChar(Buffer[21])[0]) )
     return(SUCCESS);

  if( Buffer[16] != 'A' && Buffer[16] != ' ' && Buffer[16] != '1' ) 
    return(SUCCESS);

  if( *First_ATOM ) {
    for( CC=0; CC<MAX_CHAIN; CC++ ) 
      strcpy(LastRes[CC],"XXXX");
    *First_ATOM = NO;
  }
  
  /* Find if chain exists */
  for (CC = 0; CC < *ChainNumber && Chain[CC] != NULL && Chain[CC]->Id[0] != Buffer[21]; CC++);

  if (CC == *ChainNumber) { 
    InitChain(&Chain[CC]); 
    Chain[CC]->Id[0] = Buffer[21];
    Chain[CC]->Id[1] = '\0';

    (*ChainNumber)++;
    }

//   char chainBuf[2] = { Buffer[21], '\0' };  /* 1 char + null */

//   for (CC = 0; CC < *ChainNumber; CC++) {
//     // Print CC
//     printf("ChainNumber=%d\n", CC);
//     if (Chain[CC] == NULL) {
//         printf("ERROR: Chain[%d] is NULL\n", CC);
//         break;
//     }

//     if (Chain[CC]->Id == NULL) {
//         printf("ERROR: Chain[%d]->Id is NULL\n", CC);
//         break;
//     }

//     printf("Comparing Chain[%d]->Id='%s' with chainBuf='%s'\n", CC, Chain[CC]->Id, chainBuf);


//     if (strcmp(Chain[CC]->Id, chainBuf) == 0) {
//         break;  /* We found the chain that matches Buffer[21] */
//     }
//   }

//   if (CC == *ChainNumber) {
//     InitChain(&Chain[CC]);

//     Chain[CC]->Id = (char *)malloc(2); /* big enough for 1 char + null */
//     if (!Chain[CC]->Id) {
//         fprintf(stderr, "Out of memory allocating chain ID\n");
//         return FAILURE; 
//     }
//     Chain[CC]->Id[0] = chainBuf[0];
//     Chain[CC]->Id[1] = chainBuf[1];

//     (*ChainNumber)++;
//   }
  else
  if( Chain[CC]->Ter == 1 ) 
    return(SUCCESS);

  if( Buffer[34] != '.' || Buffer[42] != '.' || Buffer[50] != '.' )
    return(escape(FAILURE,"File %s has no coordinates\n",Cmd->InputFile));

  
  if( Cmd->Stringent && Buffer[63] != '.')
    return(escape(FAILURE,"File %s has no temperature factor\n",Cmd->InputFile));


  SplitString(Buffer+22,Field,1);
  if( strcmp(Field[0],LastRes[CC]) ) {
    if( strcmp(LastRes[CC],"XXXX") && !FindAtom(Chain[CC],Chain[CC]->NRes,"CA",&NA) ) {
      free(Chain[CC]->Rsd[Chain[CC]->NRes]);
      Chain[CC]->NRes--;
    }
    if( strcmp(LastRes[CC],"XXXX") ) Chain[CC]->NRes++;
    NR = Chain[CC]->NRes;
    strcpy(LastRes[CC],Field[0]);
    Chain[CC]->Rsd[NR] = (RESIDUE *)ckalloc(sizeof(RESIDUE));
    strcpy(Chain[CC]->Rsd[NR]->PDB_ResNumb,LastRes[CC]);
    Chain[CC]->Rsd[NR]->NAtom = 0;
    SplitString(Buffer+17,Field,1);
    strcpy(Chain[CC]->Rsd[NR]->ResType,Field[0]);
  }
  else 
    NR = Chain[CC]->NRes;
  
  NA = Chain[CC]->Rsd[NR]->NAtom;

  if( Buffer[16] != ' ' ) {
    strcpy(Tmp,Buffer);
    Tmp[16] = ' ';
    SplitString(Tmp+12,Field,1);
  }
  else
    SplitString(Buffer+12,Field,1);
  
  r = Chain[CC]->Rsd[NR];
  strcpy(r->AtomType[NA],Field[0]);


  strcpy(Tmp,Buffer);
  Buffer[38] = ' ';
  SplitString(Tmp+30,Field,1);
  r->Coord[NA][0] = atof(Field[0]);

  strcpy(Tmp,Buffer);
  Buffer[46] = ' ';
  SplitString(Tmp+38,Field,1);
  r->Coord[NA][1] = atof(Field[0]);

  strcpy(Tmp,Buffer);
  Buffer[54] = ' ';
  SplitString(Tmp+46,Field,1);
  r->Coord[NA][2] = atof(Field[0]);

  if( Buffer[57] == '.' ) {
    strcpy(Tmp,Buffer);
    Tmp[60] = ' ';
    SplitString(Tmp+54,Field,1);
    r->Occupancy[NA] = atof(Field[0]);
  }
  else 
    r->Occupancy[NA] = -1.00;
  
  SplitString(Buffer+63,Field,1);
  r->TempFactor[NA] = atof(Field[0]);

  r->NAtom++;

  if( r->NAtom > MAX_AT_IN_RES-1 )
    return(escape(FAILURE,"File %s has too many atoms in residue %s %s %c\n",
		  Cmd->InputFile,r->ResType,r->PDB_ResNumb,Chain[CC]->Id));

  return(SUCCESS);
}
