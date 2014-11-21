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
  unsigned char   ta_lut[256];
      
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

  WORD = atoi(argv [3]);     //**************  Tama√±o de palabra    *******
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
		p_l_sss->nom_seq[aux]= 0;
	        printf("\n First sequence =  %s ", p_l_sss->nom_seq);
		sss->p_sss=p_l_sss;
   	     }
	     else {//---> add in head of list
                sss->n_seq++;
		for (i=0; i<(aux); i++)  p_l_sss->nom_seq[i]= ss2[i];
		p_l_sss->nom_seq[aux]= 0;
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
    

    /********************************************************************
    /* TODO QUITADO: HAY QUE METERLO EN FUNCION!!!!!

    //----->> Comienza ANOTACI”N EN LISTAS LIGADAS POR PALABRAS:
    puntero=0;
    for(i=0; i<(WORD-2); i++){ //---> WORD-1  primeros nucleotidos:
      leido = ta_lut[ss1[i]];
      puntero = puntero << 2;
      // puntero = puntero & MASCARA;
      puntero = puntero | leido;     
    }
          
    for(i=WORD-1; i<len1; i++){ //---> resto de nucleotidos:
      leido = ta_lut[ss1[i]];
      puntero = puntero << 2;
      puntero = puntero & MASCARA;
      puntero = puntero | leido;     
          
      // Nueva mejora: si no es el ˙ltimo introducido (inserciÛn por cabeza),
      // no hace falta recorrer toda la lista,
      // puesto que se asegura que ya no va a estar
      if(p_lista[puntero] == NULL || p_lista[puntero]->n_secuencia != secuencia) {
	ppunt = (struct L_NODO *) malloc(sizeof(struct L_NODO));              
	if(ppunt == NULL)  {
	  printf("\n Memory allocation failed.\n");
	  exit(EXIT_FAILURE);
	}
	ppunt->n_secuencia = secuencia;
	ppunt->mult = 1;
	
	ppunt->l_pos = (struct L_POS *) malloc(sizeof(struct L_POS));              
	if(ppunt->l_pos == NULL)  {
	  printf("\n Memory allocation failed.\n");
	  exit(EXIT_FAILURE);
	}
	//--> anoto posiciÛn de la palabra dentro de la secuencia:   
	ppunt->l_pos->posicion = i;
	ppunt->l_pos->sig_pos  = NULL;
	
	
	ppunt->sig = p_lista[puntero];
	p_lista[puntero] = ppunt;
	
	tabla[puntero]++; //--> numero de secuencias por palabra
      }
      else {
	++p_lista[puntero]->mult;
	//--> amplio la lista de posiciones dentro de una misma secuencia
	//--> para una palabra dada.
	
	p_pos = (struct L_POS *) malloc(sizeof(struct L_POS));              
	if(p_pos == NULL)  {
	  printf("\n Memory allocation failed.\n");
	  exit(EXIT_FAILURE);
	}
	p_pos->posicion= i;
	//--> anyado en cabeza de lista:
	p_pos->sig_pos = p_lista[puntero]->l_pos;
	p_lista[puntero]->l_pos = p_pos;
	//--> de momento, no incrementamos este valor
	// tabla[puntero]++; //--> numero de secuencias por palabra
	
      }
       
      long1++;
    }//-----> del FOR 
          
    ****************************************************************************/




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

  printf("\n");
  //  printf("\nSecuencias leidas =  %d \n", secuencia);

  printf("\nSecuencias leidas =  %d \n", sss->n_seq); 

  // printf("\nMaxima longitud de comentario =  %d \n", max_comen);  

  //======================================================================//
  //====    update the sequence data in the global strucutre WC     ======


  // Reservo la estructura de datos global:  
  int k; 
  k= WORD;
  wc = wc_new(k); 
  wc->k = k;

  // pointer to first sequence data:
  p_l_sss= sss->p_sss;

  for (i=0;i<sss->n_seq; i++){
    if(p_l_sss){
      //--> update sequence:
      wc_update(p_l_sss->nom_seq, p_l_sss->sequence, wc, ta_lut, i);
      p_l_sss= p_l_sss->next;
    }
  }






  printf ("\n==================================================");





 
  wc_display(wc);

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
