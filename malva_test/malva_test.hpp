/**
 * MALVA-TEST is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MALVA-TEST is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MALVA-TEST; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#ifndef _MALVA_TEST_HPP_
#define _MALVA_TEST_HPP_

#include <iostream>
#include <stdexcept>
#include <string.h>
#include "htslib/vcf.h"

typedef struct
{
   htsFile *bcf;
   bcf_hdr_t *header;
   bcf1_t *record;
} VCFt;

// VALUE of TOLERANCE RANGE MATCH of GQ
extern int tolerance; // DEFAULT 0

/** INPUT -> TWO FILEs VCF
 * OUTPUT -> Print of % Precision match
 **/
void compare_vcf(const char *sample_vcf, const char *geno_vcf);

/** Initialize the vcf file and create a VCFf struct for easier reading.
 * INPUT -> VCF file path
 * OUTPUT -> struct {bcf, header, record}
 **/
VCFt read_vcf(const char *vcf);

/** Compare #REF and #ALT
 * INPUT -> Two records
 * OUTPUT -> 0 (if #REF and #ALT are equal), 1 (else)
 **/
int equal_allele(bcf1_t *gr, bcf1_t *sr);

/** Compare GENOTYPES_VALUE(GT) in #DONOR (GT:GQ)
 * INPUT -> Size of GENO GT, GENO and SAMPLE VCFt
 * OUTPUT -> 0 (if GT are equal), 1 (else)
 **/
int equal_gt(const uint8_t size, VCFt geno, VCFt sample);

/** Compare GENOTYPES_QUALITY(GQ) in #DONOR (GT:GQ)
 * INPUT -> GENO and SAMPLE VCFt
 * OUTPUT -> 0 (if GQ are equal, considering the tolerance), 1 (else)
 **/
int equal_gq(VCFt geno, VCFt sample);

/** Print GENOTYPE_VALUE(GT) and QUALITY(GQ) in #DONOR (GT:GQ)
 * INPUT -> Size of GENO GT, GENO VCFt
 **/
void print_genotypes(const uint8_t size, VCFt geno);

/** Print all #ALT
 * INPUT -> Geno records Pointer
 **/
void print_alt(bcf1_t *gr);

#endif //_MALVA_TEST_HPP_