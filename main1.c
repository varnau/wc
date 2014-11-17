//==============================================================================
// MAIN1.C: 
//==============================================================================
// 
//==============================================================================


#include "wc.h"


//--> longitud string READ
#define SSTT    30000
#define SSTTMINI   80

#define LETRAS             4  // Es el numero de NUCLEOTIDOS
#define MAX_SECU       100    // NUMERO DE SECUENCIAS --> pasar a lista ligada
#define BARRA              50
#define UNKBYTE          1024
#define UNMEGA        1048576


void  veo_syntax(char *);

//----------------------------------------------------------

void main(int argc, char *argv[]) {
  unsigned char   ta_lut[256];
      
  FILE      *f_in, *f_out;
  char      caracter;
  int       i, j, aux;
 
  wc_t   *wc;
  struct l_node *ppunt, *paux;
  struct l_pos  *p_pos;
  
  double    aux_double;
  int       palabras1=0;
  long int  aux_long;
  int       WORD;         // = k
  long int  MEMORIA;      // 4  exp k
  int       MASCARA=0;    // Mascara 
  int       bases, long1;

  int register   puntero;

  int register   leido;
  long int       nbytes1;

  int            pos;
  char           ficha[82];
 
  struct l_numss  *num_secu;
  
  int       secuencia=0, max_comen=0;
 
  FILE      *f_sal;

 //----------------------------------------------------------
  char  ch, *ss1, *ss2, *saux;  
  int   len1, len2, llaux;
  int   n_secu;
 
  int   AA[SSTT]; //--> array auxiliar.
  int   **ta_co;  //--> tabla de coincidencias entre READs.
  int   n_tc;     //--> contador maximo del arrar ta_co[]
 
  //  int   PALABRA;
 
  //--> para al final quitar
  char  linea[SSTT];


  //--------------------------------------------------------
  //---------------------------------------------------------
  if (argc < 4) {
    printf("\nSYNTAX: \n");
    veo_syntax(argv[0]);
    // fflush(stdin); getchar();
    exit(1);
  }



  //==========================================================================
  //==============   Inicializaciones  iniciales  ============================

  for (i=0; i<256; i++) ta_lut[i] = 4;   // inicializo tabla lut   *********/
  ta_lut['A'] = 0;  ta_lut['a'] = 0;
  ta_lut['C'] = 1;  ta_lut['c'] = 1;
  ta_lut['G'] = 2;  ta_lut['g'] = 2;
  ta_lut['T'] = 3;  ta_lut['t'] = 3;
  // ta_lut['N'] = ENE;  ta_lut['n'] = ENE;  // Nuevo en 4-IV-2004
  // ta_lut['>'] = MAYOR; 
  
  //===================================================================
  //============    FICHERO de entrada   ==============================

  strcpy(ficha, argv[1]);  //------ FICHERO entrada  1  : FILE-1.FASTQ
 
  if ( (f_in=fopen(ficha,"r"))!=NULL)  //-->>  Pruebo con el nombre tal cual
    printf ("\n File_in : %s", ficha);
  else {
    strcpy(ficha, argv[1]);
    strcat(ficha, ".fastq");                  //  o Fichero.FASTQ
    if ( (f_in=fopen(ficha,"r"))==NULL) {
      printf ("\n Don't find  the File-1 %s.\n", argv[1]);
      exit(0);
    }
    else
      printf ("\n File_in : %s", ficha);
  }

  fseek(f_in, 0, SEEK_END);       //------- NOU: per saber el tamany del fitxer
  nbytes1 = ftell(f_in);
  printf("\t  --> Size = %ld bytes\n", nbytes1);
  fseek(f_in, 0, SEEK_SET);

  //===================================================================
  //============    FICHERO de salida    ==============================

  strcpy(ficha, argv[2]);  //------ FICHERO entrada  1  : FILE-1.FASTQ
  f_out=fopen(ficha,"w");
  printf ("\n File_OUT: %s", ficha);

  //============    KK : tamanyo de palabra   ==========================

  WORD = atoi(argv [3]);     //**************  Tamaño de palabra    *******
  if ( (WORD<6) || (WORD>15) ) {
    printf ("\n K = [6, 15]\n");
    exit(0);
  }

  //========================================================================
  //**********  Reservo memoria para las READ 
  ss1 =  (char *) malloc(SSTT * sizeof(char));
  if (ss1== NULL){
    printf("\n memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }
        
  //**********  Reservo memoria para nombre de secuencias
  ss2 =  (char *) malloc(SSTTMINI * sizeof(char));
  if (ss2== NULL){
    printf("\n memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }        


  int k;
  k = WORD;

  // Reservo la estructura de datos global:
  wc = wc_new(k);



  printf ("\n==================================================");











  wc_display(wc);

  wc_free(wc);
}
//----------------------------------------------------------
//----------------------------------------------------------



void  veo_syntax(char *ss){
  printf("\n==========================================================\n");
  printf("%s  FASTQ_File  Out_File  K ", ss);
  printf("\n-----------\n  k       : Word Size.");
  printf("\n==========================================================\n");

 return;
}
