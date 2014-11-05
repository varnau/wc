#ifndef _WC_H
#define _WC_H

//----------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//----------------------------------------------------------
// structs
//----------------------------------------------------------
// l_node l_node_t
typedef struct l_node {
  struct l_node *next;
  int    n_seq;
  short  mult; 
  struct l_pos *l_pos;
} l_node_t;
  
typedef struct l_pos {
  struct l_pos *next;
  int pos;
} l_pos_t;  

typedef struct wc {
  int k;
  int num_words;
  l_node_t **nodes;
} wc_t;

//----------------------------------------------------------

typedef struct wc_cmp {
} wc_cmp_t;

//----------------------------------------------------------

typedef struct wc_table {
} wc_table_t;

//----------------------------------------------------------
// wc functions
//----------------------------------------------------------

wc_t *wc_new(int k);
void wc_free(wc_t *wc);

void wc_display(wc_t *wc);
void wc_upate(char *id, char *seq, wc_t *wc);
int wc_compare(char *id1, char *id2, wc_t *wc);
wc_cmp_t *wc_full_compare(char *id1, char *id2, wc_t *wc);

//----------------------------------------------------------
// wc_cmp functions
//----------------------------------------------------------

//----------------------------------------------------------
// wc_table functions
//----------------------------------------------------------

//----------------------------------------------------------
#endif // end of _WC_H
