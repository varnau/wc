#include "wc.h"


//----------------------------------------------------------
// wc functions
//----------------------------------------------------------


// ---------------------------------------------------------
// wc_new: create a data structure
// -------
wc_t *wc_new(int k) {
  int size;
  wc_t *wc = (wc_t *) calloc(1, sizeof(wc_t));
  if (wc) {
    size = pow(4, k);
    wc->nodes = calloc(size, sizeof(l_node_t));


    if (wc->nodes) {
      wc->k = k;
      wc->num_words = size;
    } else {
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
      free(frec);
    }
    free(wc);
  }
}

//----------------------------------------------------------

void wc_display(wc_t *wc) {
  printf("K value = %i\n", wc->k);
  printf("Num. words = %i\n", wc->num_words);
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
