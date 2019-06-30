#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"

#define VERSION "0.0.1"

int main(int argc,char** argv) {
    fprintf(stderr, "Usage is: %s inputBam prob1 outputBam1 prob2 outputBam2 (" VERSION ")\n", argv[0]);
    fprintf(stderr, "  prob1 and prob2 are the probabilities that the reads in inputBam are in outputBam1 and outputBam2, respectively.\n", argv[0]);
    
    samFile *in = NULL;
    bam1_t *b= NULL;
    bam_hdr_t *header = NULL;
    samFile *out1 = NULL;
    samFile *out2 = NULL; 
    
    if (argc < 5) return -1;
    in = sam_open(argv[1], "r");
    if (in==NULL) return -1;
    if ((header = sam_hdr_read(in)) == 0) return -1;
    b = bam_init1();
    
    double prob1 = atof(argv[2]);
    double prob2 = atof(argv[4]);
    out1 = sam_open(argv[3], "wb2");
    out2 = sam_open(argv[5], "wb2");
    
    int sam_write_ret = 0;
    sam_write_ret = sam_hdr_write(out1, header);
    assert(sam_write_ret >= 0);
    sam_write_ret = sam_hdr_write(out2, header);
    assert(sam_write_ret >= 0);
    
    // the largest Mersenne prime that has exact binary representation in the double primitive type.
    // uint32_t hash_max = (0x1UL << 31UL) - 1UL; 
    // fprintf(stderr, "hash_max = %u\n", hash_max);
    while (sam_read1(in, header, b) >= 0) {
        const uint32_t subsam_seed = 0;
        /* https://github.com/samtools/samtools/blob/52bd699e0a75b6d39504d1f1bdb38dacdc903921/sam_view.c#L93
         * uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ settings->subsam_seed);
         * if ((double)(k&0xffffff) / 0x1000000 >= settings->subsam_frac) return 1;
         * */
        uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ subsam_seed);
        double randfrac = (double)(k&0xffffff) / 0x1000000;
        if (randfrac < prob1) {
            sam_write_ret = sam_write1(out1, header, b);
            assert(sam_write_ret);
        }
        if (randfrac >= 1.0 - prob2) {
            sam_write_ret = sam_write1(out2, header, b);
            assert(sam_write_ret);
        }
        /*
        uint32_t strhash = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ 0) % hash_max;
        if ((0.5 + (double)strhash) / ((double)hash_max) < prob1) {
            sam_write_ret = sam_write1(out1, header, b);
            assert(sam_write_ret);
        }
        if ((0.5 + (double)strhash) / ((double)hash_max) > 1.0 - prob2) {
            sam_write_ret = sam_write1(out2, header, b);
            assert(sam_write_ret);
        }
        */
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out1);
    sam_close(out2);
    return 0;
}

