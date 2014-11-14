#ifndef _WC_H
#define _WC_H

//----------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//---------------------------------------------------------

/******************************************
// 4 nucleotides
#define LETRAS  4
#define BARRA              50
#define UNKBYTE          1024
#define UNMEGA        1048576                                                                  

// max size of strings
#define SSTT    80000  // for a read
#define STMINI     98  // head of a read
*********************************************/

// To read FASTQ
#define MAYOR         62
#define CR            10
#define ARROBA        64
#define MAS           43
#define AND           38
#define TAB            9
   
//----------------------------------------------------------
// structs : l_node_t   l_pos_t  wc_t  l_ss_t
//----------------------------------------------------------

// ------------------------------------------------------------------
// l_node_t: nodos con las secuencias que contienen una palabra dada
//-------
// next  : siguiente secuencia
// n_sec : sequence number
// mult  : numero de veces que una palabra esta en una secuencia dada
// l_pos : punter a las posiciones de una palabra en una secuencia
// -----------------------------------------------------------------
typedef struct l_node {
  struct l_node *next;
  int    n_seq;
  short  mult; 
  struct l_pos *l_pos;
} l_node_t;
//------------------------------------------------------------------  


// ----------------------------------------------------------------
// l_pos_t : lista de posiciones de una palabra en una secuencia.
// --------
// next  : next position
// pos   : position
// 
typedef struct l_pos {
  struct l_pos *next;
  int pos;
} l_pos_t;  
// ----------------------------------------------------------------


// ----------------------------------------------------------------
// wc_t : global structure of data
// ---------
// k    : size of word
// num_word : 4 exp(k)= number of words
// nodes: pointer to nodes
// frec : number of sequences for a word
// ---------
typedef struct wc {
  int k;
  int num_words;
  l_node_t **nodes;
  size_t   *table;
} wc_t;
// ----------------------------------------------------------------

// ----------------------------------------------------------------
// l_ss_t : sequence information
// ------
// nom_seq   : first lines of sequence (comment)
// sequence  : suquence of bases
// size_seq  : sequence size
// offset    : sequence position in the FASTQ file
// ------
typedef struct l_numss { 
  char *nom_seq;
  char *sequence; //--> NOT: too much information
  size_t  size_seq;
  size_t  offset;    //--> not used
} l_ss_t;
// ---------------------------------------------------------------


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
int  wc_compare(char *id1, char *id2, wc_t *wc);
wc_cmp_t *wc_full_compare(char *id1, char *id2, wc_t *wc);

//----------------------------------------------------------
// wc_cmp functions
//----------------------------------------------------------


//----------------------------------------------------------
// wc_table functions
//----------------------------------------------------------


//----------------------------------------------------------
#endif // end of _WC_H
