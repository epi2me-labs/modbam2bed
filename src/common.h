#ifndef _MODBAMBED_COMMON_H
#define _MODBAMBED_COMMON_H

#include <stdint.h>


typedef struct mod_base {
    char *name;
    char *abbrev;
    char base;
    char code;
} mod_base;

static const size_t n_mod_bases = 10;
static const mod_base mod_bases[] = {
    {"5-methylcytosine", "5mC", 'C', 'm'},
    {"5-hydroxymethylcytosine", "5hmC", 'C', 'h'},
    {"5-formylcytosine", "5fC", 'C', 'f'},
    {"5-carboxylcytosine", "5caC", 'C', 'c'},
    {"5-hydroxymethyluracil", "5hmU", 'T', 'g'},
    {"5-formyluracil", "5fU", 'T', 'e'},
    {"5-carboxyluracil", "5caU", 'T', 'b'},
    {"6-methyladenine", "6mA", 'A', 'a'},
    {"8-Oxoguanine", "8oxoG", 'G', 'o'},
    {"Xanthosine", "Xao", 'T', 'n'},
};
static const mod_base default_mod_base = {"5-methylcytosine", "5mC", 'C', 'm'};


/** Simple integer min/max
 * @param a
 * @param b
 *
 * @returns the min/max of a and b
 *
 */
static inline int max ( int a, int b ) { return a > b ? a : b; }
static inline int min ( int a, int b ) { return a < b ? a : b; }


/** Allocates zero-initialised memory with a message on failure.
 *
 *  @param num number of elements to allocate.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xalloc(size_t num, size_t size, char* msg);


/** Retrieves a substring.
 *
 *  @param string input string.
 *  @param postion start position of substring.
 *  @param length length of substring required.
 *  @returns string pointer.
 *
 */
char *substring(char *string, int position, int length);

#endif
