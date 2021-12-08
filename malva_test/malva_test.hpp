/**
 * MALVA - genotyping by Mapping-free ALternate-allele detection of known VAriants
 * Copyright (C) 2019  Giulia Bernardini, Luca Denti, Marco Previtali
 *
 * This file is part of MALVA.
 *
 * MALVA is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MALVA is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MALVA; see the file LICENSE. If not, see
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

#endif //_MALVA_TEST_HPP_