#ifndef SUMMARY_METADATA_H_
#define SUMMARY_METADATA_H_

#include <string>

#include <stdio.h>
#include <stdint.h>

#include "khash.h"
#include "utils.h"

// The class summarizes the overall mapping metadata 

namespace chromap {

enum SummaryMetadataField {
  SUMMARY_METADATA_TOTAL = 0,
  SUMMARY_METADATA_DUP,
  SUMMARY_METADATA_MAPPED,
  SUMMARY_METADATA_LOWMAPQ,
	SUMMARY_METADATA_CACHEHIT,
  SUMMARY_METADATA_CARDINALITY,
  //SUMMARY_METADATA_CACHEHIT_FWD, 
  SUMMARY_METADATA_FIELDS
};

struct _barcodeSummaryMetadata {
  int counts[SUMMARY_METADATA_FIELDS];
  _barcodeSummaryMetadata() {
    memset(counts, 0, sizeof(int) * SUMMARY_METADATA_FIELDS);
  }
};


KHASH_MAP_INIT_INT64(k64_barcode_metadata, struct _barcodeSummaryMetadata)

class SummaryMetadata {
 public:
  SummaryMetadata() {
    barcode_metadata_ = kh_init(k64_barcode_metadata);
    barcode_length_ = 16;
  }
  ~SummaryMetadata() {
    kh_destroy(k64_barcode_metadata, barcode_metadata_);
  }

  void Output(const char *filename, std::vector<double> score_coeffs, bool output_peak_info) {
    FILE *fp = fopen(filename, "w");

    // change header based on command line options
    if (!output_peak_info)
      fprintf(fp, "barcode,total,duplicate,unmapped,lowmapq,cachehit,fric,score\n");
    else  
      fprintf(fp, "barcode,total,duplicate,unmapped,lowmapq,cachehit,numslots,fric,score\n");
    
    // iterate through all barcodes
    khiter_t k;
    for (k = kh_begin(barcode_metadata_); k != kh_end(barcode_metadata_); ++k) {
      if (kh_exist(barcode_metadata_, k)) 
      {
        // define variables to store values
        double fric = 0.0;
        double num_dup = kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_DUP]; 
        double num_unmapped = kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_TOTAL] - kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_MAPPED];
        double num_lowmapq = kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_LOWMAPQ];

        // print barcode as string
        fprintf(fp, "%s", Seed2Sequence(kh_key(barcode_metadata_, k), barcode_length_).c_str());
        int i;

        // prints out total # of reads, # duplicate, # unmapped, # lowmapq, # cache hit, and cardinality 
        for (i = 0; i < SUMMARY_METADATA_FIELDS; ++i) 
        {
          if (i == SUMMARY_METADATA_CARDINALITY && !output_peak_info) {continue;}
          if (i != SUMMARY_METADATA_MAPPED) {
            fprintf(fp, ",%d", kh_value(barcode_metadata_, k).counts[i]);
          } else { 
            // the print is for unmapped
            fprintf(fp, ",%d", kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_TOTAL]
                - kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_MAPPED]);
          }
        }
        // print out fric: fraction in cache
        fric = (double)kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_CACHEHIT] / (double)kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_MAPPED];
        fprintf(fp, ",%.5lf", fric);

        // compute the chromap score and print it out
        double chromap_score = score_coeffs[0] + /* constant */
                               (score_coeffs[1] * fric) +
                               (score_coeffs[2] * num_dup) +
                               (score_coeffs[3] * num_unmapped)  +
                               (score_coeffs[4] * num_lowmapq);
        fprintf(fp, ",%.5lf", chromap_score);
        fprintf(fp, "\n");
      }
    }
    fclose(fp);
  }

  void UpdateCount(uint64_t barcode, int type, int change) {
    int khash_return_code;
    khiter_t barcode_metadata_iter = kh_put(k64_barcode_metadata, barcode_metadata_, barcode, &khash_return_code);
    if (khash_return_code) {
      struct _barcodeSummaryMetadata nb;
      kh_value(barcode_metadata_, barcode_metadata_iter) = nb;
    }
    kh_value(barcode_metadata_, barcode_metadata_iter).counts[type] += change;
  }

  void SetBarcodeLength(int l) {
    barcode_length_ = l;
  }

  // In SAM format for paired-end data, some count will be counted twice
  void AdjustPairedEndOverCount() {
    khiter_t k;
    for (k = kh_begin(barcode_metadata_); k != kh_end(barcode_metadata_); ++k)
      if (kh_exist(barcode_metadata_, k)) {
        kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_DUP] /= 2 ;
        kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_LOWMAPQ] /= 2 ;
        kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_MAPPED] /= 2 ;
      } 
  }

 private:
  khash_t(k64_barcode_metadata) *barcode_metadata_;    
  int barcode_length_;

  std::string Seed2Sequence(uint64_t seed, uint32_t seed_length) const {
    std::string sequence;
    sequence.reserve(seed_length);
    uint64_t mask_ = 3;
    for (uint32_t i = 0; i < seed_length; ++i) {
      sequence.push_back(
          Uint8ToChar((seed >> ((seed_length - 1 - i) * 2)) & mask_));
    }
    return sequence;
  }
};

} // namespace chromap

#endif
