//
// Created by Marco on 24/05/2021.
//

#pragma once

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

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <utility>

#include <unordered_map>
#include <set>
#include <string>
#include <vector>

#include <math.h>
#include <zlib.h>

#include "htslib/hts_log.h"
#include "kmc/kmc_file.h"
#include "kseq.h"
#include "htslib/vcf.h"

#include "argument_parser.hpp"
#include "bloom_filter.hpp"
#include "var_block.hpp"
#include "kmap.hpp"

using namespace std;

static const char* MALVA_IDX_SUFFIX = ".malvax";

void pelapsed(const string &s = "", const bool rollback = false);

/**
 * Method to add kmers to the bloom filter
 **/
void add_kmers_to_bf(BF &bf, KMAP &ref_bf, const VK_GROUP &kmers);

/**
 * Method to compute and store the coverages of the alleles of the
 * variants of a var_block. It uses the coverages stored in the bloom
 * filters/map.
 **/
void set_coverages(BF &bf, KMAP &ref_bf, VB &vb, const VK_GROUP &kmers/*, const float &cap*/);

/**
 * Method to clean and print VCF header. It adds GT and GQ FORMAT,
 * removes all samples, and adds donor sample.
 **/
void print_cleaned_header(bcf_hdr_t *vcf_header, const bool verbose);

// ---------------------------------------------------------------------------

int index_main(int argc, char *argv[]);
int call_main(int argc, char *argv[]);
