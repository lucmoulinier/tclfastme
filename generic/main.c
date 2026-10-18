#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <tcl.h>
#include "main.h"
#include "graph.h"
#include "newick.h"

#include "interface_utilities.h"
#include "interface_options.h"

#ifndef DLLEXPORT
#define DLLEXPORT __declspec(dllexport)
#endif 

/* Tcl 8.6 compatibility for the Tcl 9 64-bit size API */
#ifndef TCL_SIZE_MAX
# define TCL_SIZE_MAX INT_MAX
# ifndef Tcl_Size
    typedef int Tcl_Size;
# endif
# define TCL_SIZE_MODIFIER ""
#endif


/*functions from newick.c*/
tree *loadNewickTree(FILE *ifile, int numLeaves);
Tcl_Obj *NewickPrintTree(tree *T);
tree *detrifurcate(tree *T);
void partitionSizes(tree *T);
double GetTreeLength (tree *T);

/*functions from bme.c*/
tree *BMEaddSpecies(tree *T,node *v, double **D, double **A); 
void makeBMEAveragesTable(tree *T, double **D, double **A);
void assignBMEWeights(tree *T, double **A);

/*from gme.c*/
tree *GMEaddSpecies(tree *T,node *v, double **D, double **A);
void makeOLSAveragesTable(tree *T, double **D, double **A);
void assignOLSWeights(tree *T, double **A);
void assignAllSizeFields(tree *T);

/*from BIONJ.c*/
tree *bionj (double **D, set *S, int n, boolean isNJ);

/*from inputs.c*/
double **loadMatrix(FILE *ifile, int *size, set *S);
void CopyMatrix (double **D, float **M, int n);
//void compareSets(tree *T, set *S, FILE *ofile);
void compareSets(tree *T, set *S);
void freeMatrix(double **D, int size);
void ReplaceForbiddenChar (char *str, char tofind, char replace);
void ReplaceForbiddenChars (char *str, char replace);

/*from NNI.c*/
//void NNI(tree *T, double **avgDistArray, int *count);
void NNI(tree *T, double **avgDistArray, int *count, double **D, int numSpecies);

/*from bNNI.c*/
//void bNNI(tree *T, double **avgDistArray, int *count);
void bNNI(tree *T, double **avgDistArray, int *count, double **D, int numSpecies);

/*from graph.c*/
void freeSet(set *S);
void freeTree(tree *T);
node *copyNode(node *v);

/* Pour Tcl/Tk init */
int Fastme_ObjCmd(ClientData clientData,
		  Tcl_Interp *interp,
		  int objc,
		  Tcl_Obj *const objv[]);

/* ================ */
/* = Programme ==== */
/* ================ */

void chooseSettings(int argc, char **argv, int *btype, int *ntype, int *wtype, int *numDataSets, char **filenames)
{
  int counter = 1;
  int i;

  while (counter < argc)
    {
      switch(argv[counter][1])
	{	
	case 't':
	  counter++;
	  if (NULL != argv[counter])
	    strncpy(filenames[2], argv[counter], MAX_FILE_NAME_LENGTH);
	  else
	    {
	      fprintf(stderr,"Error: -t flag requires argument.\n");
	      exit(EXIT_FAILURE);
	    }
	  counter++;
	  *btype = USER;
	  break;
	case 'w':
	  counter++;
	  switch(argv[counter][0])
	    {
	    case 'b':
	    case 'B':
	      *wtype = BAL;
	      break;
	      //	    case 'g':
	      //	    case 'G':
	    case 'o':
	    case 'O':
	      *wtype = OLS;
	      break;
	    case 'n':
	    case 'N':
	      *wtype = NONE;
	      break;
	    default:
	      fprintf(stderr,"Unknown argument to -w option: please");
	      fprintf(stderr," use (b)alanced, (O)LS or (n)one\n");
	      exit(EXIT_FAILURE);
	    }
	  counter++;
	  break;
	case 'b':
	  counter++;
	  switch(argv[counter][0])
	    {
	    case 'b':
	    case 'B':
	      if (argv[counter][1] == 'i' || argv[counter][1] == 'I')
	        *btype = BIONJ;
	      else if (argv[counter][1] == 'm' || argv[counter][1] == 'M')
	        *btype = BAL;
	      break;
	      //	    case 'g':
	      //case 'G':
	    case 'o':
	    case 'O':
	      *btype = OLS;
	      break;
	    case 'n':
	    case 'N':
	      *btype = NJ;
	      break;

	    case 'i':
	    case 'I':
	      *btype = BIONJ;
	      break;
	    default:
	      fprintf(stderr,"Unknown argument to -b option: please");
	      fprintf(stderr," use (b)alanced_GME, (O)LS_GME, B(I)ONJ or (N)J\n");
	      exit(EXIT_FAILURE);
	    }
	  counter++;
	  break;
	case 'n':
	  counter++;
	  *numDataSets = i = 0;
	  while (argv[counter][i])
	    *numDataSets = 10* (*numDataSets) + (argv[counter][i++] - '0');
	  counter++;
	  break;
	case 's':
	  counter++;
	  switch(argv[counter][0])
	    {
	    case 'b':
	    case 'B':
	      *ntype = BAL;
	      break;
	      //	    case 'g':
	      //	    case 'G':
	    case 'o':
	    case 'O':
	      *ntype = OLS;
	      break;
	    case 'n':
	    case 'N':
	      *ntype = NONE;
	      break;
	    default:
	      fprintf(stderr,"Unknown argument to -s option: please");
	      fprintf(stderr," use (b)alanced_NNI, (O)LS_NNI, or none\n");
	      exit(EXIT_FAILURE);
	    }
	  counter++;
	  break;
	  /*	case 'v':
		verbose = TRUE;
		counter++;
		break;*/
	case 'h':
	default:
	  Usage();
	  /*	  fprintf(stderr,"Usage: fastme -binostv\n");
		  fprintf(stderr,"-b specify method for building initial tree: ");
		  fprintf(stderr,"GME or BME(default).\n");
		  fprintf(stderr,"-i filename of distance matrix\n");
		  fprintf(stderr,"-n number of trees/matrices input (default = 1)\n");
		  fprintf(stderr,"-o filename for tree output\n");
		  fprintf(stderr,"-s specify type of tree swapping (NNIs): ");
		  fprintf(stderr,"(b)alanced, (O)LS, or (n)one. (Default is balanced.)\n");
		  fprintf(stderr,"-t (optional) filename of starting tree topology\n");

		  fprintf(stderr,"-v for verbose output\n");
		  fprintf(stderr,"-w (b)alanced or (O)LS branch lengths (if not doing NNIs on input topology) \n");
		  fprintf(stderr,"-help to get this message\n");
		  exit(0);
	  */	}
    }
  if ( (*btype == OLS || *btype == BAL) && *wtype == NONE)
    Usage();
}


double **initDoubleMatrix(int d)
{
  int i,j;
  double **A;
  A = (double **) malloc(d*sizeof(double *));
  for(i=0;i<d;i++)
    {
      A[i] = (double *) malloc(d*sizeof(double));
      for(j=0;j<d;j++)
	A[i][j] = 0.0;
    }
  return(A);
}


void printOptions (int ntype, int numSets, int btype, int wtype, FILE *f, char *open_mode)
{
  if (open_mode[0] == 'w') {
    fprintf (f, "\n - FastME %s - \n\n\n", VERSION);
    fprintf (f, "\nPapers to be cited:\n");
    fprintf (f, "\nFastME algorithms (balanced and OLS GME, balanced and OLS NNI, balanced and OLS branch length optimization):");
    fprintf (f, "\n\tDesper R., Gascuel O. 2002. Fast and accurate phylogeny reconstruction algorithms based on the minimum-evolution principle.");
    fprintf (f, "\n\tJournal of Computational Biology. 9(5):687-705.");
    fprintf (f, "\nBIONJ algorithm:");
    fprintf (f, "\n\tGascuel O. 1997. BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data.");
    fprintf (f, "\n\tMolecular Biology and Evolution, 14(7):685-695");
    fprintf (f, "\nNJ algorithm:");
    fprintf (f, "\n\tSaitou N., Nei M. 1987. The neighbor-joining method: a new method for reconstructing phylogenetic trees.");
    fprintf (f, "\n\tMolecular Biology and Evolution, 1987 4(4):406-25\n\n");
  }

  fprintf (f, "-------------------------------------------------------------------------------\n"); 
  fprintf (f, "Settings for this run:\n\n");

  switch (btype) {
  case OLS:
    fprintf (f, "  B "
	     " Initial tree: build method (balanced_GME, OLS_GME, BIONJ or NJ) or user tree "
	     " %-15s \n", "OLS_GME");
    break;
  case BAL: 
    fprintf (f, "  B "
	     " Initial tree: build method (balanced_GME, OLS_GME, BIONJ or NJ) or user tree "
	     " %-15s \n", "balanced_GME");
    break;
  case NJ:
    fprintf (f, "  B "
	     " Initial tree: build method (balanced_GME, OLS_GME, BIONJ or NJ) or user tree "
	     " %-15s \n", "NJ");
    break;
  case BIONJ:
    fprintf (f, "  B "
	     " Initial tree: build method (balanced_GME, OLS_GME, BIONJ or NJ) or user tree "
	     " %-15s \n", "BIONJ");
    break;
  case USER:
    fprintf (f, "  B "
	     " Initial tree: build method (balanced_GME, OLS_GME, BIONJ or NJ) or user tree "
	     " %-15s \n", "USER");
    break;
  default:
    break;
  }
  //----------------------------------------------------------------------//
  fprintf (f, "  N "
	   "                                                     Number of matrices/trees "
	   " %-15d \n", numSets);
  //----------------------------------------------------------------------//
  switch (ntype) {
  case OLS:
    fprintf (f, "  S "
	     "                          Type of tree swapping (NNI) (balanced, OLS or none) "
	     " %-15s \n", "OLS");
    break;
  case BAL:
    fprintf (f, "  S "
	     "                          Type of tree swapping (NNI) (balanced, OLS or none) "
	     " %-15s \n", "balanced");
    break;
  case NONE:
    fprintf (f, "  S "
	     "                          Type of tree swapping (NNI) (balanced, OLS or none) "
	     " %-15s \n", "none");
    //----------------------------------------------------------------------//
    switch (wtype) {
    case OLS:
      fprintf (f, "  W "
	       "                                      Branch lengths assigned to the topology "
	       " %-15s \n", "OLS");
      break;
    case BAL:
      fprintf (f, "  W "
	       "                                      Branch lengths assigned to the topology "
	       " %-15s \n", "balanced");
      break;
    case NONE:
      fprintf (f, "  W "
	       "                                      Branch lengths assigned to the topology "
	       " %-15s \n", "none");
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  fprintf (f, "\n");
  return;
}

/*
  Tcl interface for FastME. 
  This routine is a mix between the 'main' routine and the 'chooseSettings' routine.
*/
int DLLEXPORT Tclfastme_Init(Tcl_Interp *interp)
{
  if (Tcl_InitStubs(interp, "8.6-", 0) == NULL) {
    return TCL_ERROR;
  }
  if (Tcl_PkgProvide(interp, PACKAGE_NAME, PACKAGE_VERSION) != TCL_OK) {
    return TCL_ERROR;
  }
  
  // Creation de la commande Tcl
  Tcl_CreateObjCommand(interp, "fastme", Fastme_ObjCmd, NULL, (Tcl_CmdDeleteProc*)NULL);
  
  return TCL_OK;
}


/*
  Parse the value of -init_tree, -branch_lengths or -swap.
  Values are matched on their first letters, case-insensitively.
  The Tcl string is only read, never modified.
*/
static int GetMethodFromObj(Tcl_Interp *interp, Tcl_Obj *objPtr,
			    const char *option, int allowNJ, int allowNone,
			    int *typePtr)
{
  const char *buf = Tcl_GetString(objPtr);

  if (strncasecmp(buf,"bal",3) == 0) {
    *typePtr = BAL;
  } else if (strncasecmp(buf,"ols",3) == 0) {
    *typePtr = OLS;
  } else if (allowNJ && strncasecmp(buf,"bio",3) == 0) {
    *typePtr = BIONJ;
  } else if (allowNJ && strncasecmp(buf,"nj",2) == 0) {
    *typePtr = NJ;
  } else if (allowNone && strncasecmp(buf,"non",3) == 0) {
    *typePtr = NONE;
  } else {
    Tcl_SetObjResult(interp, Tcl_ObjPrintf(
      "bad value \"%s\" for %s: must be %s", buf, option,
      allowNJ ? "bal, ols, bionj or nj" : "bal, ols or none"));
    return TCL_ERROR;
  }
  return TCL_OK;
}

/*
  The FastME code keeps its state in globals (nodeCount, edgeCount).
  Serialize calls so that interpreters in different threads are safe.
*/
TCL_DECLARE_MUTEX(fastmeMutex)

int Fastme_ObjCmd(ClientData clientData,
		  Tcl_Interp *interp,
		  int objc,
		  Tcl_Obj *const objv[])
{
  int i, j;
  double **D = NULL, **A = NULL;
  Tcl_Obj **nameObjv, **distObjv;
  Tcl_Size nameObjc, distObjc;
  node *v;
  set *species, *slooper;
  tree *T;
  const char *nom;
  char noeud[NODE_LABEL_LENGTH];
  node *addNode;
  int numSpecies;
  int index;
  int nniCount = 0;
  int ntype, btype, wtype;

  static const char *const optionsTcl[] = {
    "-init_tree", "-branch_lengths", "-swap", NULL};
  enum optionsTcl {
    FME_INIT_TREE,FME_BRANCH,FME_SWAP
  };

  if (objc < 3 || ((objc - 3) % 2) != 0) {
    Tcl_WrongNumArgs(interp, 1, objv,
      "names distances ?-init_tree method? ?-branch_lengths method? ?-swap method?");
    return TCL_ERROR;
  }

  /* Defaults: balanced GME tree, balanced NNI */
  ntype = BAL;
  wtype = btype = BAL;

  for (i = 3; i < objc; i += 2) {
    if (Tcl_GetIndexFromObj(interp, objv[i], optionsTcl, "option", 0, &index) != TCL_OK) {
      return TCL_ERROR;
    }
    switch ((enum optionsTcl) index) {
    case FME_BRANCH:
      if (GetMethodFromObj(interp, objv[i+1], "-branch_lengths", 0, 1, &wtype) != TCL_OK)
	return TCL_ERROR;
      break;
    case FME_INIT_TREE:
      if (GetMethodFromObj(interp, objv[i+1], "-init_tree", 1, 0, &btype) != TCL_OK)
	return TCL_ERROR;
      break;
    case FME_SWAP:
      if (GetMethodFromObj(interp, objv[i+1], "-swap", 0, 1, &ntype) != TCL_OK)
	return TCL_ERROR;
      break;
    }
  }

  if (Tcl_ListObjGetElements(interp, objv[1], &nameObjc, &nameObjv) != TCL_OK
      || Tcl_ListObjGetElements(interp, objv[2], &distObjc, &distObjv) != TCL_OK) {
    return TCL_ERROR;
  }
  if (nameObjc == 0) {
    return TCL_OK;
  }
  /* FastME indexes species with int; n*n doubles must also be addressable */
  if (nameObjc > 65535) {
    Tcl_SetObjResult(interp, Tcl_ObjPrintf(
      "too many species (%" TCL_SIZE_MODIFIER "d)", nameObjc));
    return TCL_ERROR;
  }
  numSpecies = (int) nameObjc;
  if ((Tcl_WideInt) distObjc != (Tcl_WideInt) numSpecies * numSpecies) {
    Tcl_SetObjResult(interp, Tcl_ObjPrintf(
      "distance list has %" TCL_LL_MODIFIER "d elements, expected %d x %d = %"
      TCL_LL_MODIFIER "d", (Tcl_WideInt) distObjc, numSpecies, numSpecies,
      (Tcl_WideInt) numSpecies * numSpecies));
    return TCL_ERROR;
  }

  /* Read the distances before allocating the tree structures */
  D = initDoubleMatrix(numSpecies);
  for (i = 0; i < numSpecies; i++) {
    for (j = 0; j < numSpecies; j++) {
      if (Tcl_GetDoubleFromObj(interp, distObjv[(Tcl_Size)i*numSpecies+j], &D[i][j]) != TCL_OK) {
	freeMatrix(D,numSpecies);
	return TCL_ERROR;
      }
    }
  }

  /* Species names: truncated to fit a node label, and made Newick-safe */
  species = (set *) malloc(sizeof(set));
  species->firstNode = NULL;
  species->secondNode = NULL;
  for (i = 0; i < numSpecies; i++) {
    nom = Tcl_GetString(nameObjv[i]);
    strncpy(noeud,nom,NODE_LABEL_LENGTH-1);
    noeud[NODE_LABEL_LENGTH-1] = '\0';
    ReplaceForbiddenChars (noeud, '_');
    for (slooper = species; slooper != NULL && slooper->firstNode != NULL;
	 slooper = slooper->secondNode) {
      if (strcmp(slooper->firstNode->label, noeud) == 0) {
	Tcl_SetObjResult(interp, Tcl_ObjPrintf(
	  "duplicate species name \"%s\"", noeud));
	freeSet(species);
	freeMatrix(D,numSpecies);
	return TCL_ERROR;
      }
    }
    v = makeNewNode(noeud,-1);
    v->index2 = i;
    species = addToSet(v,species);
  }

  /* A single species: the tree is just the leaf */
  if (numSpecies == 1) {
    Tcl_SetObjResult(interp, Tcl_ObjPrintf("%s;\n", species->firstNode->label));
    freeSet(species);
    freeMatrix(D,numSpecies);
    return TCL_OK;
  }

  /* BIONJ/NJ need at least 3 taxa; with fewer the tree is unique anyway */
  if (numSpecies < 3 && (btype == BIONJ || btype == NJ)) {
    btype = BAL;
  }

  Tcl_MutexLock(&fastmeMutex);

  A = initDoubleMatrix(2*numSpecies-2);
  T = NULL;

  switch(btype) {
  case OLS:
    for(slooper = species; NULL != slooper; slooper = slooper->secondNode) {
      addNode = copyNode(slooper->firstNode);
      T = GMEaddSpecies(T,addNode,D,A);
    }
    break;
  case BAL:
    for(slooper = species; NULL != slooper; slooper = slooper->secondNode) {
      addNode = copyNode(slooper->firstNode);
      T = BMEaddSpecies(T,addNode,D,A);
    }
    break;
  case NJ:
  case BIONJ:
    T = bionj (D, species, numSpecies, (btype == NJ));
    T = detrifurcate(T);
    compareSets(T,species);
    partitionSizes(T);
    break;
  }

  switch(ntype) {
  case OLS:
    if (OLS != btype)
      assignAllSizeFields(T);
    makeOLSAveragesTable(T,D,A);
    assignOLSWeights(T,A);
    NNI(T,A,&nniCount,D,numSpecies);
    assignOLSWeights(T,A);
    break;
  case BAL:
    if (BAL != btype)
      makeBMEAveragesTable(T,D,A);
    assignBMEWeights(T,A);
    bNNI(T,A,&nniCount,D,numSpecies);
    assignBMEWeights(T,A);
    break;
  case NONE:
    switch(wtype) {
    case OLS:
      if (OLS != btype)
	assignAllSizeFields(T);
      makeOLSAveragesTable(T,D,A);
      assignOLSWeights(T,A);
      break;
    case BAL:
      if (BAL != btype)
	makeBMEAveragesTable(T,D,A);
      assignBMEWeights(T,A);
      break;
    }
    break;
  }

  Tcl_SetObjResult(interp, NewickPrintTree(T));

  freeTree(T);
  freeMatrix(A,2*numSpecies - 2);
  Tcl_MutexUnlock(&fastmeMutex);

  freeMatrix(D,numSpecies);
  freeSet(species);

  return TCL_OK;
}


/*
int tralala_main(int argc, char **argv)
{
  int i;
  int setCounter = 0;
  FILE *ifile1, *ifile2, *ofile, *statfile, *bionjfile;
  double **D, **A;
  set *species, *slooper;
  node *addNode;
  char **filenames;
  tree *T;
  int numSpecies;
  int numSets = 1;
  int count = 0;
  int nniCount = 0;
  int ntype, btype, wtype;
  option *input;
  char *open_mode = "a";
  struct stat tmpbuf;

  ntype = BAL;
  wtype = btype = BAL;
  //  verbose = FALSE;
  filenames = (char **) malloc (4*sizeof(char *));
  T = NULL;
  for(i=0;i<4;i++)
    filenames[i] = (char *) malloc(MAX_FILE_NAME_LENGTH*sizeof(char));
    
  if (argc == 1) {	// PHYLIP like interface
    fflush (NULL);
    input = (option *) Get_Input (argc, argv);
    filenames[0] = input->mat_file;
    filenames[1] = input->tree_file;
    filenames[2] = input->topo_file;
    filenames[3] = input->stat_file;
    if (input->tree_file_open_mode == 0)
      open_mode = "r";
    else if (input->tree_file_open_mode == 1)
      open_mode = "w";
    else if (input->tree_file_open_mode == 2)
      open_mode = "a";
      
    numSets = input->nb_matrix;
    if (input->branch_length == 'b')
      wtype = BAL;
    else if (input->branch_length == 'O')
      wtype = OLS;
    else if (input->branch_length == 'n')
      wtype = NONE;
    if (input->NNI == 'b')
      ntype = BAL;
    else if (input->NNI == 'O')
      ntype = OLS;
    else if (input->NNI == 'n')
      ntype = NONE;
    if (strncmp (input->method, "balanced", 8) == 0)
      btype = BAL;
    else if (strncmp (input->method, "OLS", 3) == 0)
      btype = OLS;
    else if (strncmp (input->method, "NJ", 2) == 0)
      btype = NJ;
    else if (strncmp (input->method, "BIONJ", 5) == 0)
      btype = BIONJ;
    if (input->use_topo_file == TRUE)
      btype = USER;
  }
  else {			// Command line interface
    chooseSettings(argc,argv,&btype,&ntype,&wtype,&numSets,filenames);
    // if output file doesn't exist, open it in writing mode
    if (stat (filenames[1], &tmpbuf) == -1) {
      open_mode = "w";
    }
  }

  ifile1 = fopen(filenames[0],"r");
  if (!ifile1)
    {
      fprintf(stderr,"Error opening input file %s\n",filenames[0]);
      exit(EXIT_FAILURE);
    }
  bionjfile = fopen(filenames[0],"r");
  if (!bionjfile)
    {
      fprintf(stderr,"Error opening input file %s\n",filenames[0]);
      exit(EXIT_FAILURE);
    }
  ofile = fopen(filenames[1],open_mode);
  if (!ofile)
    {
      fprintf(stderr,"Error opening output file %s\n",filenames[1]);
      exit(EXIT_FAILURE);
    }
  statfile = fopen(filenames[3],open_mode);
  if (!statfile)
    {
      fprintf(stderr,"Error opening output statistic file %s\n",filenames[3]);
      exit(EXIT_FAILURE);
    }
  printOptions (ntype, numSets, btype, wtype, statfile, open_mode);

  fprintf(stdout,"numSpecies %d\n",numSpecies);
  fflush(stdout);
  if (USER == btype)
    ifile2 = fopen(filenames[2],"r");
    // lulu
  while ( setCounter < numSets ) {
    setCounter++;
    species = (set *) malloc(sizeof(set));
    species->firstNode = NULL;
    species->secondNode = NULL;

    D = loadMatrix(ifile1,&numSpecies,species);
    fprintf(stdout,"numSpecies %d\n",numSpecies);
    fflush(stdout);
    A = initDoubleMatrix(2*numSpecies-2);
    fprintf (statfile,"Dataset n° %d\n", setCounter);
    switch(btype)
      {
      case USER:
	if (!ifile2)
	  {
	    fprintf(stderr,"Error opening input file %s\n",filenames[2]);
	    exit(EXIT_FAILURE);
	  }
	T = loadNewickTree(ifile2,numSpecies);
	T = detrifurcate(T);
	//compareSets(T,species,ofile);
	compareSets(T,species);
	partitionSizes(T);
	fclose (ifile2);
	break;
      case OLS:
	for(slooper = species; NULL != slooper; slooper = slooper->secondNode)
	  {
	    addNode = copyNode(slooper->firstNode);
	    T = GMEaddSpecies(T,addNode,D,A);
	  }
	break;
      case BAL:
	for(slooper = species; NULL != slooper; slooper = slooper->secondNode)
	  {
	    addNode = copyNode(slooper->firstNode);
	    T = BMEaddSpecies(T,addNode,D,A);
	  }
	break;
      case NJ:
        T = bionj (bionjfile, true);
        T = detrifurcate(T);
        //compareSets(T,species,ofile);
        compareSets(T,species);
	partitionSizes(T);
        break;
      case BIONJ:
        T = bionj (bionjfile, false);
        T = detrifurcate(T);
        //compareSets(T,species,ofile);
        compareSets(T,species);
	partitionSizes(T);
        break;
      }

    switch(ntype)
      {
      case OLS:
	if (OLS != btype)
	  assignAllSizeFields(T);
	makeOLSAveragesTable(T,D,A);
	// assign branch length to T to compute tree length for statistics
	assignOLSWeights(T,A);
	fprintf (statfile,"\tBefore performing NNI(s), tree length is: %f\n", GetTreeLength (T));
	//NNI(T,A,&nniCount);
	NNI(T,A,&nniCount,D,numSpecies);
	fprintf (statfile, "\tNumber of NNI(s) performed: %d\n", nniCount);
	assignOLSWeights(T,A);
	break;
      case BAL:
	if (BAL != btype)
	  makeBMEAveragesTable(T,D,A);
	// assign branch length to T to compute tree length for statistics
	assignBMEWeights(T,A);
	fprintf (statfile,"\tBefore performing NNI(s), tree length is: %f\n", GetTreeLength (T));
	//bNNI(T,A,&nniCount);
	bNNI(T,A,&nniCount,D,numSpecies);
	fprintf (statfile, "\tNumber of NNI(s) performed: %d\n", nniCount);
	assignBMEWeights(T,A);
	break;
      case NONE:
	switch(wtype)
	  {
	  case OLS:
	    if (OLS != btype)
	      assignAllSizeFields(T);
	    makeOLSAveragesTable(T,D,A);
	    assignOLSWeights(T,A);
	    break;
	  case BAL:
	    if (BAL != btype)
	      makeBMEAveragesTable(T,D,A);
	    assignBMEWeights(T,A);
	    break;
	  case NONE:
	    break;
	  default:
	    fprintf(stderr,"Error in program: variable 'btype' has illegal ");
	    fprintf(stderr,"value %d.\n",btype);
	    exit(EXIT_FAILURE);
	  }
	break;
      default:
	fprintf(stderr,"Error in program: variable 'ntype' has illegal ");
	fprintf(stderr,"value %d.\n",ntype);
	exit(EXIT_FAILURE);
      }

    fprintf (statfile, "\tTree length is: %f\n\n", GetTreeLength (T));

    NewickPrintTree(T,ofile);
    freeMatrix(D,numSpecies);
    freeMatrix(A,2*numSpecies - 2);
    freeSet(species);
    freeTree(T);
    T = NULL;
    //    setCounter++;
    //    if ((verbose) && (ntype))
    //      printf("Performed %d NNIs on data set %d\n",nniCount,setCounter);
    nniCount = 0;
  }
  if (argc == 1)
    Free_Input (input);
  fclose (ofile);
  fclose (statfile);
  fclose (bionjfile);
  fclose (ifile1);
  exit (0);
}
*/


