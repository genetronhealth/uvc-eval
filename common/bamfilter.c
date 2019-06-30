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

int main(int argc,char** argv) {
    
    samFile *in = NULL;
    bam1_t *b= NULL;
    bam_hdr_t *header = NULL;
    samFile *out1 = NULL;
    samFile *out2 = NULL; 
    
    if (argc < 4) {
        fprintf(stderr, "Usage is: %s input-bam pass-output-bam fail-output-bam <pass-condition>\n", argv[0]);
        fprintf(stderr, "  pass-condition is nonzero-length cigar operation length by default.\n", argv[0]);
        return -1;
    }
    in = sam_open(argv[1], "r");
    if (NULL == in) {
        fprintf(stderr, "Input error for file %s at initialization.\n", argv[1]);
        return -1;
    }
    
    if ((header = sam_hdr_read(in)) == 0) {
        fprintf(stderr, "Input error for file %s at header.\n", argv[1]);
        return -1;
    }
    b = bam_init1();
    
    out1 = sam_open(argv[2], "wbu");
    out2 = sam_open(argv[3], "wbu");
    
    int sam_write_ret = 0;
    sam_write_ret = sam_hdr_write(out1, header);
    assert(sam_write_ret >= 0);
    sam_write_ret = sam_hdr_write(out2, header);
    assert(sam_write_ret >= 0);
    
    int64_t nreads = 0;
    int64_t nreads_failed = 0;
    int64_t nreads_thres = (100 * 1000);
    while (sam_read1(in, header, b) >= 0) {
        nreads++;
        const uint32_t *cigar = bam_get_cigar(b);
        int is_pass = 1;
        for (int i = 0; i < b->core.n_cigar; i++) {
            if (0 == bam_cigar_oplen(cigar[i])) { is_pass = 0; }
        }
        if (is_pass) {
            sam_write_ret = sam_write1(out1, header, b);
            assert(sam_write_ret || !fprintf(stderr, "Writting of passed read %s at tid %d pos %d failed!\n", bam_get_qname(b), b->core.tid, b->core.pos));
        } else {
            nreads_failed++;
            sam_write_ret = sam_write1(out2, header, b);
            assert(sam_write_ret || !fprintf(stderr, "Writting of failed read %s at tid %d pos %d failed!\n", bam_get_qname(b), b->core.tid, b->core.pos));
        }
        if (nreads > nreads_thres) {
            fprintf(stderr, "Processed %ld reads and %ld processed reads failed. Processed tid %d pos %d\n", nreads, nreads_failed, b->core.tid, b->core.pos);
            nreads_thres = nreads_thres * 11 / 10 + (100*1000);
        }
    }
    
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out1);
    sam_close(out2);
    return 0;
}

