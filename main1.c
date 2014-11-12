#include <stdio.h>
#include <stdlib.h>

#include "wc.h"

void  veo_syntax(char *);

//----------------------------------------------------------

void main(int arg, char *argv[]) {
  unsigned char   ta_lut[256];
      
  FILE      *f_in, *f_out;
  char      caracter;
  int       i, j, aux;
 
  int  *tabla;  // number of sequence
  struct l_node  **, *ppunt, *paux;
  struct l_pos   *p_pos;
  
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
 
  struct l_ss_t  *num_secu;
  
  int       secuencia=0, max_comen=0;
 
  FILE      *f_sal;

 //----------------------------------------------------------
  char  ch, *ss1, *ss2, *saux;  
  int   len1, len2, llaux;
  int   n_secu;
 
  int   AA[SSTT]; //--> array auxiliar.
  int   **ta_co;  //--> tabla de coincidencias entre READs.
  int   n_tc;     //--> contador maximo del arrar ta_co[]
 
  int   PALABRA;
 
  //--> para al final quitar
  char  linea[SSTT];


  //--------------------------------------------------------
  //---------------------------------------------------------
  if (argc < 4) {
    printf("\nSYNTAX: \n\n");
    veo_syntax(argv[0]);
    fflush(stdin); getchar();
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

  PALABRA = atoi(argv [3]);     //**************  Tamaño de palabra    *******
  if ( (PALABRA<6) || (PALABRA>15) ) {
    printf ("\n KK = [6, 15]\n");
    exit(0);
  }






  int k = 10;



  wc_t *wc = wc_new(k);






  wc_display(wc);

  wc_free(wc);
}
//----------------------------------------------------------
//----------------------------------------------------------
