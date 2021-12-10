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

int equal_allele(bcf1_t *gr, bcf1_t *sr){
    for(int i = 0; i < (gr->n_allele); ++i){
        if(strcmp(gr->d.allele[i], sr->d.allele[i]) != 0){
            return 1;
        }
    }
    return 0;
}

int equal_gt(const uint8_t size, VCFt geno, VCFt sample){
    //EXTRACT xEX -> X|A 
    uint8_t geno_n1 = bcf_get_fmt(geno.header,geno.record,"GT")->p[0];
    uint8_t sample_n1 = bcf_get_fmt(sample.header,sample.record,"GT")->p[0];
    //HAPLOID
    if((size == 1) && (geno_n1 == sample_n1)){
        return 0;
    }
    //NORMAL then EXTRACT yEY -> B|Y
    if((size == 2) && (geno_n1 == sample_n1)){
        if(bcf_get_fmt(geno.header,geno.record,"GT")->p[1]
           == bcf_get_fmt(sample.header,sample.record,"GT")->p[1]){
               return 0;
           }
    }
    return 1;
}

int equal_gq(VCFt geno, VCFt sample){
    //EXTRACT GQ GENO
    int geno_gq = (int)bcf_get_fmt(geno.header,geno.record,"GQ")->p[0];
    //EXTRACT GQ SAMPLE
    int sample_gq = (int)bcf_get_fmt(sample.header,sample.record,"GQ")->p[0];
    //CALCULATE THE TOLERANCE RANGE [min, max]
    int min, max;
    if( (geno_gq-tolerance) < 0 ){
        min = 0;
    } else{
        min = geno_gq - tolerance;
    }
    if( (geno_gq+tolerance) > 100 ){
        max = 100;
    } else{
        max = geno_gq + tolerance;
    }
    //CHECK IF SAMPLE_GQ is in the TOLERANCE RANGE
    if( (min<=sample_gq) && (sample_gq<= max) ){
        return 0;
    }
    return 1;
}

void print_genotypes(const uint8_t size, VCFt geno){
    std::cerr << " #DONOR ";
    //PRINT GT
    uint8_t n1 = bcf_get_fmt(geno.header,geno.record,"GT")->p[0];
    if(n1 == 2){ // 0|Y
        std::cerr << "0";
    }
    if(n1 == 4){ // 1|Y
        std::cerr << "1";
    }
    if(size == 2){
        uint8_t n2 = bcf_get_fmt(geno.header,geno.record,"GT")->p[1];
        if(n2 == 2){ // X|0
            std::cerr << "|0";
        }
        if(n2 == 4){ // X|1
            std::cerr << "|1";
        }
    }
    //PRINT GQ
    std::cerr << ":" << unsigned(bcf_get_fmt(geno.header,geno.record,"GQ")->p[0]);
}

void print_alt(bcf1_t *gr){
    if((gr->n_allele) > 1){ //if have ALTs
        std::cerr << " #ALT ";
    }
    int i = 1;
    while(i < (gr->n_allele)){
        std::cerr << gr->d.allele[i];
        ++i;
        if(i < (gr->n_allele)){
            std::cerr << ",";
        }
    }
}

void compare_vcf(const char* geno_vcf, const char* sample_vcf){
    //PRINT USED FILEs
    std::cerr << std::endl << "Compare \"" << geno_vcf << "\" with \"" << sample_vcf << "\"" << std::endl << std::endl;
    
    //LOAD GENO
    VCFt geno = read_vcf(geno_vcf);
    //LOAD SAMPLE
    VCFt sample = read_vcf(sample_vcf);
    
    double compari = 0; //number of comparisons
    double match = 0; //number of records matched
    
    /***
        * 1. read GENO_RECORD and SAMPLE_RECORD
        * 2. check GENO_RECORD are covered in the SAMPLE_RECORD
        * (same: #CHROM #POS #ID #REF #ALT)
        * 3. compari++
        * 4. if(record covered) match++ else (print it)
    ***/
    while( (bcf_read(geno.bcf, geno.header, geno.record) == 0) 
            && (bcf_read(sample.bcf, sample.header, sample.record) == 0) ) {
        compari++;
        
        //unpack GENO_RECORD for read REF,ALT,INFO,etc 
        bcf_unpack(geno.record,BCF_UN_ALL);

        //unpack SAMPLE_RECORD for read REF,ALT,INFO,etc 
        bcf_unpack(sample.record,BCF_UN_ALL);

        //DEBUG: print all records comparisons (#CHROM #POS #REF)
        if(DEBUG){
            std::cerr << "Record Geno #CHROM " << geno.record->rid+1 << ", #POS " << geno.record->pos+1 << ", #REF " << geno.record->d.allele[0]<< std::endl;            
            std::cerr << "Record Sample #CHROM " << sample.record->rid+1 << ", #POS " << sample.record->pos+1 << ", #REF " << sample.record->d.allele[0]<< std::endl << std::endl;
        }

        /***
            *1) COMPARE #CHROM (int64_t)
            *2) COMPARE #POS (int64_t)
            *3) COMPARE #ID (char* str)
            *4) COMPARE #REF and #ALT (char* str)
            *5) COMPARE GENOTYPES_VALUE_SIZE (uint8_t) in #DONOR (example: size 2 -> 0|1 ; size 1(haplide) -> 0) 
            *6) COMPARE GENOTYPES_VALUE(GT) in #DONOR (example: 0|1) 
            *7) COMPARE GENOTYPES_QUALITY(GQ) in #DONOR (:100)
        ***/
        if( (geno.record->rid+1 == sample.record->rid+1) 
            && (geno.record->pos+1 == sample.record->pos+1)
            && (strcmp(geno.record->d.id, sample.record->d.id) == 0)
            && (equal_allele(geno.record, sample.record) == 0) 
            && (bcf_get_fmt(geno.header, geno.record,"GT")->size == bcf_get_fmt(sample.header, sample.record,"GT")->size)
            && (equal_gt(bcf_get_fmt(geno.header,geno.record,"GT")->size, geno, sample) == 0) 
            && (equal_gq(geno, sample) == 0) ){   
            match++;
        }else{
            std::cerr << "<< NOT FOUND: " << 
            "#CHROM " << geno.record->rid+1 << 
            " #POS " << geno.record->pos+1 << 
            " #ID " << geno.record->d.id <<
            " #REF " << geno.record->d.allele[0];
            print_alt(geno.record);
            print_genotypes(bcf_get_fmt(geno.header,geno.record,"GT")->size, geno);
            std::cerr << " >>" << std::endl;
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
    std::cerr << std::endl << "Records Matched: " << match << ", Records Processed: " << compari << std::endl;
    std::cerr << "Value of Precision: " << 100*(match/compari) << "%" << std::endl << std::endl;
}