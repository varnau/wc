//==============================================================================
// MAIN1.C: 
//==============================================================================
// 
//==============================================================================

#include "wc.h"

/**************************************************
//--> longitud string READ
#define SSTT    30000
#define SSTTMINI   80

#define LETRAS             4  // Es el numero de NUCLEOTIDOS
#define MAX_SECU       100    // NUMERO DE SECUENCIAS --> pasar a lista ligada
#define BARRA              50
#define UNKBYTE          1024
#define UNMEGA        1048576
**************************************************/


void  veo_syntax(char *);

//----------------------------------------------------------

void main(int argc, char *argv[]) {
  //  unsigned char   ta_lut[256];
      
  FILE      *f_in, *f_out;
  char      caracter;
  int       i, j, aux;

  // ---------> new data structure: 
  wc_t   *wc;
  struct l_node *ppunt, *paux;
  struct l_pos  *p_pos;

  sss_t  *sss;
  l_ss_t *p_l_sss;
  // ----------------------------
  
  double    aux_double;
  int       palabras1=0;
  long int  aux_long;
  int       WORD;         // = k
  long int  MEMORIA;      // 4  exp k
  //  int       MASCARA=0;    // Mascara 
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
  ss2 =  (char *) malloc(STMINI * sizeof(char));
  if (ss2== NULL){
    printf("\n memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }        
  //**************************************
  // Reservo memoria para estructura de datos de secuencias FASTQ:
                                                                    

  //----  auxiliar  memory to strings (sequences)

  ss1 =  (char *) malloc(SSTT * sizeof(char));
  if (ss1== NULL){
    printf("\n memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }
        
  //**********  Reservo memoria para nombre de secuencias
  ss2 =  (char *) malloc(STMINI * sizeof(char));
  if (ss2== NULL){
    printf("\n memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }        

  // ---------> New data structure to sequences of FASTQ file:

  sss = sss_new();

  sss->n_seq=0;
  sss->file = argv[1];


  // **************************************************************************
  // ******  LEO FICHERO_FASTQ   y genero tabla de frecuencias   **************

  printf("\n\n------- COMIENZA EL ANALISIS  -------------- (pulsa tecla)\n"); 
  fflush(stdin); getchar();
  long1=0;                           // icincializo tamanyo FILE-1
  secuencia = 0;  //--> comenzamos con secuencia 1: se incrementa en comentario

  while (!feof(f_in)) {
    //----->> Leo primera linea:
    ch = fgetc(f_in);
    if(ch!=ARROBA){
      printf("\n\n-------> FIN DEL ANALISIS\n");  
      break;
    }          
              
    len2=0;          
    do {    
      ch = fgetc(f_in);
      ss2[len2]=ch;
      len2++;
      // printf("\n %d  : %c ", (int)ch, ch);
      // printf("%c", ch);
    }
    while(ch!=CR);   
    ss2[len2]= 0;
    
    // fflush(stdin); getchar();
 
    // ES VIEJO: QUITAR
    // strcpy(num_secu[secuencia].nom, ss2);
    // printf("\n....... %s ......", num_secu[secuencia].nom);

    // ----> incorporo datos a nueva estructura:
    aux = strlen(ss2);
    if (aux){
       p_l_sss = (l_ss_t *) calloc(1, sizeof (l_ss_t)); 
       if (p_l_sss==NULL){
           printf("\n Memory allocation failed.\n");
           exit(EXIT_FAILURE); 
       } else{
	     p_l_sss->nom_seq= calloc(1, aux+1);

             if (secuencia==0){//--> first sequence
	        sss->n_seq=1;
	        for (i=0; i<(aux); i++)  p_l_sss->nom_seq[i]= ss2[i];
		p_l_sss->nom_seq[aux-1]= 0;  //--> CR to end to sequence
		p_l_sss->num=0;               //--> number of de first sequence

	        printf("\n First sequence =  %s ", p_l_sss->nom_seq);
		sss->p_sss=p_l_sss;
   	     }
	     else {//---> add in head of list
                sss->n_seq++;
		for (i=0; i<(aux); i++)  p_l_sss->nom_seq[i]= ss2[i];
		p_l_sss->nom_seq[aux-1]= 0;  //--> CR to end to sequence
		p_l_sss->num= secuencia;      //--> number of sequence

                // strcpy( p_l_sss->nom_seq, ss2);
		printf("\n Next  sequence =  %s ", p_l_sss->nom_seq);
                p_l_sss->next=sss->p_sss;
	        sss->p_sss= p_l_sss;   
  	      }
	     //  fflush(stdin); getchar(); 
            }
    }
    else{
       printf("\n String error.\n");
       exit(EXIT_FAILURE);
    }
    
    
    // fflush(stdin); getchar();
    
    //----->> Leo SEGUNDA linea: SECUENCIA DE NUCLEOTIDOS      
    len1=0;  
    ch = fgetc(f_in); //--> leo primer nucleotido!!    
    do {
      ss1[len1]=ch;
      len1++;
      ch = fgetc(f_in);
      // printf("\n %d  : %c ", (int)ch, ch);
      // printf("%c", ch);
    }
    while(ch!=CR);      
    ss1[len1]= 0;
    
    //----> QUITAR Y REEMPLAZA:
    //    num_secu[secuencia].tam = len1;
    //  printf("\n len= %d", len1);
    
    //----->> anotamos en nueva estructura de datos:
    aux= strlen(ss1); 
    p_l_sss->sequence= calloc(1, aux+1); 
    if (p_l_sss->sequence){
        p_l_sss->size_seq= aux;    
	for (i=0; i<(aux); i++)  p_l_sss->sequence[i]= ss1[i]; 
	p_l_sss->sequence[aux] =0;
        printf("\t-->  len= %d", aux); 
        printf("\n %s", p_l_sss->sequence);
    } else {
      printf("\n Memory allocation failed in SSS.\n");
      exit(EXIT_FAILURE);
    } 
    // fflush(stdin); getchar();
    

    //----->> Leo tercera linea:
    ch = fgetc(f_in);
    if(ch!=MAS){
      printf("\n\n-------> FIN DEL ANALISIS\n");  
      break;
    }           
    do {    
      ch = fgetc(f_in);
      // printf("\n %d  : %c ", (int)ch, ch);
      // printf("%c ", ch);
    }
    while(ch!=CR);   
    // fflush(stdin); getchar();
    
    //----->> Leo cuarta linea:      
    do {    
      ch = fgetc(f_in);
      // printf("%c", ch);
      // printf("\n %d  : %c ", (int)ch, ch);
    }
    while(ch!=CR);   
    // fflush(stdin); getchar();
    secuencia++;
    
  }
 
 
  // printf("\n\n-------> FIN DEL ANALISIS!!!!\n"); 
  // fflush(stdin); getchar();
  
  fclose(f_in);    //--> close DNA_File (FASTQ file) !!!

  // printf("\nSecuencias leidas =  %d \n", secuencia);
  printf("\n Secuencias leidas =  %d \n", sss->n_seq); 

  // printf("\nMaxima longitud de comentario =  %d \n", max_comen);  

  //======================================================================
  //====    update the sequence data in the global strucutre WC     ======


  // Reservo la estructura de datos global:  
  int k; 
  k= WORD;
  wc = wc_new(k, secuencia); 
  fflush(stdin); getchar();  

  wc->k = k;
  wc->n_seq=  sss->n_seq;  

  // pointer to first sequence data:
  p_l_sss= sss->p_sss;

  for (i=0;i<sss->n_seq; i++){
    if(p_l_sss){
      //--> update sequence:
      //      wc_update(p_l_sss->nom_seq, p_l_sss->sequence, wc, i);

      //--> MAL: tenia que llamarlas al reves, por lo de la insercion en cabeza
      //  wc_update(p_l_sss->nom_seq, p_l_sss->sequence, wc, sss->n_seq-i-1);
      wc_update(p_l_sss->nom_seq, p_l_sss->sequence, wc, p_l_sss->num);

      p_l_sss= p_l_sss->next;
    }
  }


  wc_display(wc); 

  printf ("\n=================================================="); 


  //==================================================================
  //=== Genero la Tabla de palabras coincidentes entre secuencias

  printf("\n Llamo a wc_table \n");
  //  fflush(stdin); getchar();
  wc_table(wc);

  printf("\n Llamo a wc_display_table \n");  
  fflush(stdin); getchar();
  wc_display_table(wc);

  printf ("\n==================================================");

  printf("\n Llamo a wc_display_table_FILE \n");
  fflush(stdin); getchar();
  wc_display_table_FILE(wc, sss, f_out);


  printf ("\n==================================================");


  printf("\n Llamo a wc_compare_by_index \n");
  fflush(stdin); getchar();
  int  index1, index2;
  printf("\n Dame los dos indices : ");
  scanf("%d %d", &index1, &index2);


  wc_compare_by_index(index1, index2, wc);
 
  printf ("\n=================================================="); 

  printf("\n FREE  Memory ... \n");

  wc_free(wc);
  sss_free(sss);

  return;
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
