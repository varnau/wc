#include "wc.h"

//----------------------------------------------------------
// wc functions
//----------------------------------------------------------

// ---------------------------------------------------------
// wc_new: create a data structure
// -------
wc_t *wc_new(int k) {
  size_t size;
  int    i;

  wc_t *wc = (wc_t *) calloc(1, sizeof(wc_t));
  if (wc) {
    size = (size_t) pow(4, k);
    // memory for the word frecuency table
    wc->table = calloc (size, sizeof(size_t));
    // Memory for the nodes: array of pointes!! 
    wc->nodes = (l_node_t **) calloc(size, sizeof(l_node_t *));

    if (wc->nodes && wc->table) {
      wc->k = k;
      wc->num_words = size;
    } else {
      if (wc->table!=NULL) free(wc->table);
      if (wc->nodes!=NULL) free(wc->nodes); 
      free(wc);
    }
  }
  return wc;
}
//----------------------------------------------------------

// --------------------------------------------------------
// wc_free : free memory for a data structure
// --------
void wc_free(wc_t *wc) {
  if (wc) {
    l_node_t *node, *next_node;
    l_pos_t  *pos,  *next_pos;

   // printf("\n n_words = %ld ", wc->num_words);

    for (int i = 0; i < wc->num_words; i++) {
      node = wc->nodes[i];
      while(node) {
	// save next node
	next_node = node->next;

	// free node
	pos = node->l_pos; 
	while(pos) {
	  // save next pos
	  next_pos = pos->next;
	  // free pos
	  free(pos);
	  // restore pos
	  pos = next_pos;
	}
	free(node);	

	// restore node
	node = next_node;
      }
      // free frecuency words
      // free(wc->table);
    }
    // free frecuency words
    free(wc->table);
    free(wc);
  }//--> del IF(wc)
}


//-------------------------------------------------------------
// wc_update: update sequence words in global data structure WC
// ---------
void wc_update(char *id, char *seq, wc_t *wc, unsigned char *ta_lut, int secuencia) {
  int puntero;
  size_t MASCARA;
  size_t WORD;
  int leido;
  int len1, i;
 
  l_node_t *ppunt, paux;
  l_pos_t  *p_pos;


  MASCARA= (size_t) pow(4.0, wc->k); 
  MASCARA = MASCARA-1;
  printf("\n MASCARA = %ld ", MASCARA);
  WORD= wc->k;


  printf("\n %s ", id);
  printf("\n %s ", seq);
  printf("\n %d ", wc->k);

  len1=strlen(seq);

  //----->> Comienza ANOTACIÓN EN LISTAS LIGADAS POR PALABRAS:
  puntero=0;
  for(i=0; i<(WORD-2); i++){ //---> WORD-1  primeros nucleotidos:
    leido = ta_lut[seq[i]];
    puntero = puntero << 2;
    // puntero = puntero & MASCARA;
    puntero = puntero | leido;     
  }
          
  for(i=WORD-1; i<len1; i++){ //---> resto de nucleotidos:
    leido = ta_lut[seq[i]];
    puntero = puntero << 2;
    puntero = puntero & MASCARA;
    puntero = puntero | leido;     
          
    // Nueva mejora: si no es el último introducido (inserción por cabeza),
    // no hace falta recorrer toda la lista,
    // puesto que se asegura que ya no va a estar

    // BORRAR:
    // if(p_lista[puntero] == NULL || p_lista[puntero]->n_secuencia != secuencia) {


    if(wc->nodes[puntero] == NULL || wc->nodes[puntero]->n_seq != secuencia) {

      ppunt = (l_node_t *) calloc(1, sizeof(l_node_t));              
      if(ppunt == NULL)  {
	printf("\n Memory allocation failed.\n");
	exit(EXIT_FAILURE);
      }
      ppunt->n_seq = secuencia;
      ppunt->mult = 1;
      
      ppunt->l_pos = (l_pos_t *) calloc(1, sizeof(l_pos_t));              
      if(ppunt->l_pos == NULL)  {
	printf("\n Memory allocation failed.\n");
	exit(EXIT_FAILURE);
      }
      //--> anoto posición de la palabra dentro de la secuencia:   
      ppunt->l_pos->pos = i;
      ppunt->l_pos->next  = NULL;
      
      
      ppunt->next = wc->nodes[puntero];
      wc->nodes[puntero] = ppunt;
      
      wc->table[puntero]++; //--> numero de secuencias por palabra
    }
    else {
      ++wc->nodes[puntero]->mult;
      //--> amplio la lista de posiciones dentro de una misma secuencia
      //--> para una palabra dada.
      
      p_pos = (l_pos_t *) calloc(1, sizeof(l_pos_t));              
      if(p_pos == NULL)  {
	printf("\n Memory allocation failed.\n");
	exit(EXIT_FAILURE);
      }
      p_pos->pos= i;
      //--> anyado en cabeza de lista:
      p_pos->next = wc->nodes[puntero]->l_pos;
      wc->nodes[puntero]->l_pos = p_pos;
      //--> de momento, no incrementamos este valor
      // tabla[puntero]++; //--> numero de secuencias por palabra
      
    }
     
    // long1++;
  }//-----> del FOR 







}

//----------------------------------------------------------

void wc_display(wc_t *wc) {
  printf("\nK value = %i\n", wc->k);
  printf("Num. words = %ld\n", wc->num_words);

  int i;
  //  for (int i = 0; i < wc->num_words; i++) {
  for (int i = 0; i < 1024; i++) {
    if (wc->nodes[i]){
	printf("\n[%3d]--> %3d ", i, wc->nodes[i]->n_seq );
	printf("\t --> %3d ",  wc->nodes[i]->l_pos->pos);
      }
  }

}



//---------------------------------------------------------
// n_sec  functions
//---------------------------------------------------------

 l_ss_t *n_sec_new(void){
 

 }


 void n_sec_free(l_ss_t *n_sec){


 }


//----------------------------------------------------------


//----------------------------------------------------------

int wc_compare(char *id1, char *id2, wc_t *wc) {
}

//----------------------------------------------------------

wc_cmp_t *wc_full_compare(char *id1, char *id2, wc_t *wc) {
}

//----------------------------------------------------------

//----------------------------------------------------------
// sequence functions
//----------------------------------------------------------

// ---------------------------------------------------------
// sss_new: create a data structure for the FASTQ FILE
// -------
sss_t *sss_new(void) {
  sss_t *sss;

  sss = (sss_t *) calloc(1, sizeof(sss_t));
  if (sss) {
     sss->p_sss =(l_ss_t *) calloc(1, sizeof(l_ss_t));
     if (sss->p_sss) sss->n_seq=0;
     else
       free (sss); 
      }

  return sss;
}

// --------------------------------------------------------
// sss_free : free memory for sequence data structure
// --------
void sss_free(sss_t *sss) {
  l_ss_t *node, *next_node; 

  if (sss) {
    if (sss->p_sss){
       node=sss->p_sss;
       next_node = node->next;
       while(node){
	 free(node->nom_seq);
	 free(node->sequence);
	 free(node);
         node=next_node;
	 if(node)
	    next_node=node->next;
       }
    }
  }
  free(sss);

  return;
}

 
