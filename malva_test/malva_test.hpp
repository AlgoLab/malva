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

typedef struct {
    htsFile *bcf;
    bcf_hdr_t *header;
    bcf1_t *record;
} VCFt;

/* INPUT -> TWO FILEs VCF
   OUTPUT -> Print of % Precision match
 */
void compare_genotypes(const char* sample_vcf, const char* geno_vcf);

/* Initialize the vcf file and create a VCFf struct for easier reading.
   INPUT -> VCF file path
   OUTPUT -> struct {bcf, header, record}
*/
VCFt read_vcf(const char* vcf);

/* Compare #REF and #ALT
   INPUT -> Two records
   OUTPUT -> 0 (if #REF and #ALT are equal), 1 (else)
*/
int equal_allele(bcf1_t *gr, bcf1_t *sr);

#endif //_MALVA_TEST_HPP_