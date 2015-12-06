#define _GNU_SOURCE
#define  NUM_BASES 5
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define DFA_BLOCKS 2 // must be < 8

// Align state machine.
typedef struct {
   int8_t   delta;
   int8_t   gamma;
   uint16_t nexts;
} state_t;

static const int pow3[9] = {1,3,9,27,81,243,729,2187,6561};

static const int translate_convert[256] = {
   5,4,4,4,4,4,4,4,4,4, 6,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,
   4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,
   4,4,4,4,4,0,4,1,4,4, 4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,3,4,4,4,4,
   4,4,4,4,4,4,4,0,4,1, 4,4,4,2,4,4,4,4,4,4, 4,4,4,4,4,4,3,3,4,4,
   4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,
   4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,
   4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,
   4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,
   4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4
   };

int       parse          (const char *, char *);
int       compute_ident  (uint8_t *, char *, int, int);
int       partial_delta  (int, int);
state_t * compute_dfa    (int n);

int main(int argc, char *argv[])
{
   if (argc != 4) {
      fprintf(stderr, "usage: adfa <dist> <pattern> <file>\n");
      return EXIT_FAILURE;
   }

   // Distance.
   int dist = atoi(argv[1]);
   if (dist < 0) {
      fprintf(stderr, "error invalid distance\n");
      return EXIT_FAILURE;
   }

   // Open file.
   FILE * fin = fopen(argv[3], "r");
   if (fin == NULL) {
      fprintf(stderr, "error opening file %s\n", argv[2]);
      return EXIT_FAILURE;
   }

   // Parse pattern.
   char * keys = malloc(strlen(argv[2]));
   int plen = parse(argv[2], keys);
   if (plen < 1) {
      fprintf(stderr, "pattern error\n");
      return EXIT_FAILURE;
   }

   if (dist >= plen) {
      fprintf(stderr, "error distance geq than pattern length\n");
      return EXIT_FAILURE;
   }

   // Precompute identity vectors.
   int align_blocks = (plen+DFA_BLOCKS-1)/DFA_BLOCKS;
   int delta_depth  = plen % DFA_BLOCKS;
   uint8_t * ident = malloc(align_blocks*NUM_BASES);
   compute_ident(ident, keys, plen, DFA_BLOCKS);

   // Compute dfa.
   state_t * dfa = compute_dfa(DFA_BLOCKS);
   size_t s_sz = (1<<DFA_BLOCKS)*3;
   size_t m_sz = 3;

   // Align!
   size_t bsize = 100, lineno = 0;
   char * line = malloc(bsize);
   ssize_t nbytes = 0;
   uint8_t * align = malloc(align_blocks);

   while((nbytes = getline(&line, &bsize, fin)) > 0) {
      lineno++;
      // Reset alignment (Epsilon state).
      memset(align, 0x00, align_blocks);
      // Align nucleotides.
      int last_active = dist, last_score = dist+1, hit = 1;
      for (size_t i = 0; i < nbytes; i++) {
         int c = translate_convert[(int)line[i]];
         if (c < NUM_BASES) {
            uint8_t * id = ident + c*align_blocks;
            state_t s = {0,1,0};
            int active_ref = last_active+1, score = 0;
            int j = 0;
            for (; j < align_blocks; j++) {
               // Break if there's nothing interesting passed this point.
               if (score > dist && j > active_ref)
                     break;
               // Update state
               s = dfa[align[j]*s_sz + id[j]*m_sz + (s.gamma & 0x0F)];
               align[j] = s.nexts;
               // Check if this block had active sub-states.
               if (score + (s.gamma >> 4) <= dist)
                  last_active = j;
               
               // Update score.
               score += s.delta;
            }
            if (delta_depth && j == align_blocks) {
               score += partial_delta(s.nexts, delta_depth) - s.delta;
            }
            if (last_score <= dist && score > last_score && hit) {
               fprintf(stdout, "%ld:%ld,%d\n", lineno, i-1, last_score);
               hit = 0;
            } else if (score < last_score) {
               hit = 1;
            }
            last_score = score;
         } else if (last_score <= dist && hit) {
            fprintf(stdout, "%ld:%ld,%d\n", lineno, i-1, last_score);
            break;
         }
      }
   }
   return 0;
}

int
partial_delta
(
 int state,
 int depth
)
{
   int delta = 0;
   for (int i = 0; i < depth; i++) {
      int numstate = state / pow3[DFA_BLOCKS-1-i];
      delta += 1 - numstate;
      state -= numstate * pow3[DFA_BLOCKS-1-i];
   }
   return delta;
}

state_t *
compute_dfa
(
 int n
)
{
   static const state_t align_sm[3][3][2] = {{{{1,0,0},{1,0,0}},{{1,1,0},{0,0,1}},{{ 0,1,1},{-1,0,2}}},
                                             {{{1,1,0},{1,1,0}},{{1,2,0},{0,1,1}},{{ 0,2,1},{-1,1,2}}},
                                             {{{1,2,0},{1,2,0}},{{0,2,1},{0,2,1}},{{-1,2,2},{-1,2,2}}}};

   int nmatch = (1 << n);
   if (n < 1 || n > 8) return NULL;

   size_t state_width = nmatch*3;
   size_t match_width = 3;

   // Alloc DFA.
   state_t * dfa = malloc(pow3[n]*nmatch*3*sizeof(state_t));

   for (int s = 0; s < pow3[n]; s++) {
      for (int m = 0; m < nmatch; m++) {
         for (int g = 0; g < 3; g++) {
            // Compute n-step output.
            state_t st = {0,g,0};
            int numstate = s;
            int nextstate = 0;
            int delta = 0;
            int mindelta = 1;
            for (int i = 0; i < n; i++) {
               // Resolve match.
               int match = (m >> i) & 1;
               // Resolve current state.
               int state = numstate/pow3[n-1-i];
               numstate -= state*pow3[n-1-i];
               st = align_sm[state][st.gamma][match];
               nextstate += st.nexts * pow3[n-1-i];
               delta += st.delta;
               if (delta < mindelta) mindelta = delta;
            }
            dfa[s*state_width + m*match_width + g] = (state_t){delta,(st.gamma & 0x0F) | ((mindelta << 4) & 0xF0),nextstate};
         }
      }
   }
   return dfa;
}


int
parse
(
 const char * expr,
       char * keys
)
{
   // Initialize keys to 0.
   for (size_t i = 0; i < strlen(expr); i++) keys[i] = 0;

   int i = 0;
   int l = 0;
   int add = 0;
   char c, lc = 0;
   while((c = expr[i]) != 0) {
      if      (c == 'A' || c == 'a') keys[l] |= 0x01;
      else if (c == 'C' || c == 'c') keys[l] |= 0x02;
      else if (c == 'G' || c == 'g') keys[l] |= 0x04;
      else if (c == 'T' || c == 't' || c == 'U' || c == 'u') keys[l] |= 0x08;
      else if (c == 'N' || c == 'n') keys[l] |= 0x1F;
      else if (c == '[') {
         if (add) {
            return -1;
         }
         add = 1;
      }
      else if (c == ']') {
         if (!add) {
            return -1;
         }
         if (lc == '[') l--;
         add = 0;
      }
      else {
         return -1;
      }

      if (!add) l++;
      i++;
      lc = c;
   }
   
   if (add == 1) {
      return -1;
   }
   else return l;
}

int
compute_ident
(
 uint8_t * ident,
 char    * keys,
 int       wlen,
 int       blocks
 )
{
   int bytes = (wlen+blocks-1)/blocks;
   // Set memory to 'F'.
   memset(ident, 0xFF >> (8-blocks), bytes*NUM_BASES);
   for (int j = 0, o = 0; j < NUM_BASES; j++, o += bytes) {
      uint8_t base = 1 << j;
      for (int i = 0; i < wlen; i++)
         ident[o+i/blocks] &= ~(((base & keys[i]) == 0) << (i%blocks));
   }
   return 0;
}
