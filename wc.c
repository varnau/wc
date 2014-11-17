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
    l_pos_t *pos, *next_pos;
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
      free(wc->table);
    }
    free(wc);
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

void wc_display(wc_t *wc) {
  printf("K value = %i\n", wc->k);
  printf("Num. words = %ld\n", wc->num_words);
}

//----------------------------------------------------------

void wc_upate(char *id, char *seq, wc_t *wc) {
}

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
	 free(node);
	 node=next_node;
	 next_node=node->next;
       }
    }
  }
  free(sss);

  return;
}

 
