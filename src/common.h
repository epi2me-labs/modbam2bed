#ifndef _MODBAMBED_COMMON_H
#define _MODBAMBED_COMMON_H

#include <stdint.h>


typedef struct mod_base {
    char *name;
    char *abbrev;
    char base;
    int base_i; // 16bit IUPAC form A:1, C:2, G:4, T:8
    int code;   // to enable htslib ChEBI support, chars below so simplicity
} mod_base;

static const size_t n_mod_bases = 16;
static const mod_base mod_bases[] = {
    // C mods
    {"5-methylcytosine", "5mC", 'C', 2, 'm'},
    {"5-hydroxymethylcytosine", "5hmC", 'C', 2, 'h'},
    {"5-formylcytosine", "5fC", 'C', 2, 'f'},
    {"5-carboxylcytosine", "5caC", 'C', 2, 'c'},
    {"Ambiguous C modification", "modC", 'C', 2, 'C'},
    // T mods
    {"5-hydroxymethyluracil", "5hmU", 'T', 8, 'g'},
    {"5-formyluracil", "5fU", 'T', 8, 'e'},
    {"5-carboxyluracil", "5caU", 'T', 8, 'b'},
    {"Ambiguous T modification", "modT", 'T', 8, 'T'},
    // A mods
    {"6-methyladenine", "6mA", 'A', 1, 'a'},
    {"Ambiguous A modification", "modA", 'A', 1, 'A'},
    // G mods
    {"8-Oxoguanine", "8oxoG", 'G', 4, 'o'},
    {"Ambiguous G modification", "modG", 'G', 4, 'G'},
    // U mods
    {"Ambiguous U modification", "modU", 'U', 15, 'U'}, // TODO: should 15 (N) be something else?
    // N Mods
    {"Xanthosine", "Xao", 'N', 15, 'n'},
    {"Ambiguous N modification", "modN", 'N', 15, 'N'},
};
static const mod_base default_mod_base = {"5-methylcytosine", "5mC", 'C', 2, 'm'};

//0123456789ABCDEF
//=ACMGRSVTWYHKDBN  aka seq_nt16_str[]
//=TGKCYSBAWRDMHVN  comp1ement of seq_nt16_str
//084C2A6E195D3B7F
static int seqi_rc[] = { 0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15 };

static const int MAX_QUAL = 255;

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
