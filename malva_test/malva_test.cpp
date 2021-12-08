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

#include "malva_test.hpp"

const int DEBUG = 0; //=1 FOR DEBUG PRINTs

VCFt read_vcf(const char* vcf){
    VCFt out = {NULL, NULL, NULL};
        
    //VCF INIT
    htsFile *vcf_bcf = NULL;
    //OPEN FILE and SAVE FILE METADATA
    vcf_bcf = bcf_open(vcf, "r");
    if(vcf_bcf == NULL) {
        throw std::runtime_error("Unable to open vcf file.");
    }
    out.bcf = vcf_bcf;
    
    bcf_hdr_t *vcf_header = NULL;
    //EXTRACT and SAVE VCF HEALDER
    vcf_header = bcf_hdr_read(vcf_bcf);
    if(vcf_header == NULL) {
        throw std::runtime_error("Unable to read vcf header.");
    }
    out.header = vcf_header;
    
    //EXTRACT  and SAVE VCF RECORDs
    bcf1_t *vcf_record = bcf_init();
    out.record = vcf_record;
    
    return out;
}

void compare_genotypes(const char* geno_vcf, const char* sample_vcf){
    //PRINT USED FILEs
    std::cerr << "Compare \"" << geno_vcf << "\" with \"" << sample_vcf << "\"" << std::endl;
    
    //LOAD GENO
    VCFt geno = read_vcf(geno_vcf);
    //LOAD SAMPLE
    VCFt sample = read_vcf(sample_vcf);
    
    double compari = 0; //number of comparisons
    double match = 0; //number of records matched
    int first_iter = 1; //first iteration of while
    
    /***
        * 1. read GENO_RECORD and SAMPLE_RECORD
        * 2. check GENO_RECORD are covered in the SAMPLE_RECORD
        * (same: #CHROM #POS #ID #REF)
        * 3. compari++
        * 4. if(record covered) match++ else (print it)
    ***/
    while( (bcf_read(geno.bcf, geno.header, geno.record) == 0) 
            && (bcf_read(sample.bcf, sample.header, sample.record) == 0) ) {
        compari++;
        
        //if first iteration unpack the records for read REF,ALT,INFO,etc 
        if(first_iter){
            //unpack GENO_RECORD
            bcf_unpack(geno.record, BCF_UN_STR);
            bcf_unpack(geno.record, BCF_UN_INFO);

            //unpack SAMPLE_RECORD
            bcf_unpack(sample.record, BCF_UN_STR);
            bcf_unpack(sample.record, BCF_UN_INFO);
            first_iter = 0;
        }
        
        //DEBUG: print all records comparisons (#CHROM #POS #REF)
        if(DEBUG){
            std::cerr << "Record Geno #CHROM " << geno.record->rid+1 << ", #POS " << geno.record->pos+1 << ", #REF " << geno.record->d.allele[0]<< std::endl;            
            std::cerr << "Record Sample #CHROM " << sample.record->rid+1 << ", #POS " << sample.record->pos+1 << ", #REF " << sample.record->d.allele[0]<< std::endl << std::endl;
        }

        //COMPARE #CHROM (int64_t)
        //COMPARE #POS (int64_t)
        //COMPARE #ID (char* str)
        //COMPARE #REF (char* str)
        if( (geno.record->rid+1 == sample.record->rid+1) 
            && (geno.record->pos+1 == sample.record->pos+1)
            && (strcmp(geno.record->d.id, sample.record->d.id) == 0)
            && (strcmp(geno.record->d.allele[0], sample.record->d.allele[0]) == 0) ){    
            match++;
        }else{
            std::cerr << "<< RECORD NOT FOUND: " << 
            "#CHROM " << geno.record->rid+1 << 
            " #POS " << geno.record->pos+1 << 
            " #REF " << geno.record->d.allele[0] 
            << " >>" << std::endl;   
        }
    }
    
    //GENO DESTROY
    bcf_hdr_destroy(geno.header);
    bcf_destroy(geno.record); 
    bcf_close(geno.bcf);
    //SAMPLE DESTROY
    bcf_hdr_destroy(sample.header);
    bcf_destroy(sample.record); 
    bcf_close(sample.bcf);
        
    //PRINT RESULTS: Match, comparisons, % precision match of vcfs
    std::cerr << "Records Matched: " << match << " Records Processed: " << compari << std::endl;
    std::cerr << "Value of Precision: " << 100*(match/compari) << "%" << std::endl;
}