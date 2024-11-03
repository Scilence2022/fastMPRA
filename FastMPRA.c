#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <pthread.h>
#include <zlib.h>
#include <ctype.h>
#include "uthash.h"
#include "kseq.h"
#include "kthread.c"
#include "khash.h"

KHASH_SET_INIT_STR(str)

KSEQ_INIT(gzFile, gzread)

#define DEFAULT_THREADS 5
#define DEFAULT_BLOCK_SIZE 100
#define MAX_READ_LEN 1024
#define MAX_LINE 1024

// Structure for a zone
typedef struct {
    char name[64];
    int start;
    int end;
    char* left_kmer;   // Store left k-mer sequence
    char* right_kmer;  // Store right k-mer sequence
} zone_t;

// Structure for zones array
typedef struct {
    zone_t* zones;
    int count;
    int capacity;
} zones_array_t;

// Define kmer_set_t
typedef struct {
    khash_t(str) *set;
} kmer_set_t;

// Structure for design
typedef struct {
    char* seq;
    char label[64];     // Add label field
    zones_array_t zones;
    kmer_set_t kmer_set;  // Use kmer_set_t for the design
} design_t;

// Structure for a design array
typedef struct {
    design_t* designs;    // Dynamic array to store multiple designs
    int count;            // Current design count
    int capacity;         // Currently allocated capacity
} designs_array_t;

// Structure for storing reads
typedef struct {
    char *name;
    char *seq;
    char *qual;
} read_t;

// Structure for paired reads
typedef struct {
    read_t r1;
    read_t r2;
    char *assembled;  // Store assembled sequence for PE reads
} read_pair_t;

// Structure for read batch processing
typedef struct {
    read_pair_t *pairs;
    int n_pairs;
    int k;  // k-mer size
    int is_pe;  // 1 for paired-end, 0 for single-end
    int min_overlap;
    float max_mismatch_ratio;
    designs_array_t *designs;  // Add pointer to designs
} process_batch_t;

// Configuration structure
typedef struct {
    int k;
    int threads;
    int min_overlap;           // Minimum overlap length
    float max_mismatch_ratio;  // Maximum mismatch ratio

    // [Input]
    char fq1[1024];
    char fq2[1024];

    // [Output]
    char output_file[1024];    // Output file path
    int output_format;         // Output format option

    // [Advanced]
    int batch_size;            // Batch processing size
    int quality_threshold;     // Quality threshold

    // [Design]
    designs_array_t designs;  // Use designs array instead of single design
} config_t;

// Global variables for thread coordination
static int n_threads = DEFAULT_THREADS;
static int g_k = 41;  // Default k-mer size

// Function prototypes

char* extract_kmer(const char* seq, int position, int k, int direction);
void load_config(const char* filename, config_t* config);
void print_config(config_t* config);

void init_zones_array(zones_array_t* zones);
void add_zone(zones_array_t* zones, const char* name, int start, int end,
              const char* seq, int k);
void free_zones_array(zones_array_t* zones);
void init_designs_array(designs_array_t* designs);
void add_design(designs_array_t* designs, const char* seq, const char* label);
void free_designs_array(designs_array_t* designs);
void init_default_config(config_t* config);
int str_equals(const char* str1, const char* str2);
char* extract_kmer(const char* seq, int position, int k, int direction);
float compute_jaccard_index(kmer_set_t *set1, kmer_set_t *set2);
int process_zones(const char *assembled_seq, int assembled_length, int k, design_t *design);
void init_kmer_set(kmer_set_t *kmer_set);
void compute_kmer_set(const char *seq, int k, kmer_set_t *kmer_set);
int kmer_in_set(kmer_set_t *kmer_set, const char *kmer);
void free_kmer_set(kmer_set_t *kmer_set);
char* concatenate_reads(const char* seq1, const char* seq2);
int find_kmer_position(const char* seq, int seq_len, const char* kmer, int k, int start_pos);
void output_results(read_pair_t *pair, design_t *design, int is_reverse);
void init_read(read_t *r);
void free_read(read_t *r);
void reverse_complement(const char* seq, char* rev_comp_seq, int length);
void reverse_quality(const char* qual, char* rev_qual, int length);
int pe_assemble(const char* read1_seq, const char* read1_qual,
                const char* read2_seq, const char* read2_qual,
                int min_overlap, float max_mismatch_fraction,
                char* merged_seq, char* merged_qual);

// Function to trim whitespace
char* trim(char* str) {
    char* end;
    while (isspace((unsigned char)*str)) str++;
    if (*str == 0) return str;
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    return str;
}

// Function to compare strings
int str_equals(const char* str1, const char* str2) {
    return strcmp(str1, str2) == 0;
}

// Function to initialize default configuration
void init_default_config(config_t* config) {
    config->k = 41;
    config->threads = DEFAULT_THREADS;
    config->min_overlap = 12;
    config->max_mismatch_ratio = 0.1;
    config->batch_size = DEFAULT_BLOCK_SIZE;
    config->quality_threshold = 20;
    config->output_format = 0;

    // Initialize strings to empty
    config->fq1[0] = '\0';
    config->fq2[0] = '\0';
    config->output_file[0] = '\0';

    init_designs_array(&config->designs);
}

// Initialize zones array
void init_zones_array(zones_array_t* zones) {
    zones->capacity = 4;  // Initial capacity
    zones->count = 0;
    zones->zones = (zone_t*)malloc(zones->capacity * sizeof(zone_t));
}

// Add a zone to the array
void add_zone(zones_array_t* zones, const char* name, int start, int end,
              const char* seq, int k) {
    if (zones->count >= zones->capacity) {
        zones->capacity *= 2;
        zones->zones = (zone_t*)realloc(zones->zones,
                                        zones->capacity * sizeof(zone_t));
    }

    zone_t* zone = &zones->zones[zones->count];

    strncpy(zone->name, name, 63);
    zone->name[63] = '\0';
    zone->start = start;
    zone->end = end;

    // Extract left and right k-mers
    // Left k-mer from start position going k bases to the left
    zone->left_kmer = extract_kmer(seq, start, k, -1);
    // Right k-mer from end position going k bases to the right
    zone->right_kmer = extract_kmer(seq, end, k, 1);

    zones->count++;
}

// Free zones array
void free_zones_array(zones_array_t* zones) {
    for (int i = 0; i < zones->count; i++) {
        free(zones->zones[i].left_kmer);
        free(zones->zones[i].right_kmer);
    }
    free(zones->zones);
    zones->zones = NULL;
    zones->count = 0;
    zones->capacity = 0;
}

// Initialize designs array
void init_designs_array(designs_array_t* designs) {
    designs->capacity = 4;
    designs->count = 0;
    designs->designs = (design_t*)malloc(designs->capacity * sizeof(design_t));
    for (int i = 0; i < designs->capacity; i++) {
        designs->designs[i].seq = NULL;
        designs->designs[i].label[0] = '\0';  // Initialize label
        init_zones_array(&designs->designs[i].zones);
        designs->designs[i].kmer_set.set = NULL;  // Initialize kmer set
    }
}

// Add a design to the array
void add_design(designs_array_t* designs, const char* seq, const char* label) {
    if (designs->count >= designs->capacity) {
        designs->capacity *= 2;
        designs->designs = (design_t*)realloc(designs->designs,
                                              designs->capacity * sizeof(design_t));
        // Initialize newly allocated space
        for (int i = designs->count; i < designs->capacity; i++) {
            designs->designs[i].seq = NULL;
            designs->designs[i].label[0] = '\0';  // Initialize label
            init_zones_array(&designs->designs[i].zones);
            designs->designs[i].kmer_set.set = NULL;  // Initialize kmer set
        }
    }

    design_t* design = &designs->designs[designs->count];
    design->seq = strdup(seq);
    strncpy(design->label, label, 63);
    design->label[63] = '\0';

    // Initialize k-mer set
    init_kmer_set(&design->kmer_set);
    compute_kmer_set(seq, g_k, &design->kmer_set);

    designs->count++;
}

// Free designs array
void free_designs_array(designs_array_t* designs) {
    for (int i = 0; i < designs->count; i++) {
        free(designs->designs[i].seq);
        free_zones_array(&designs->designs[i].zones);
        free_kmer_set(&designs->designs[i].kmer_set);
    }
    free(designs->designs);
    designs->designs = NULL;
    designs->count = 0;
    designs->capacity = 0;
}

// Implement other required functions here...

// Function to load INI file
void load_config(const char* filename, config_t* config) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Could not open config file %s. Using defaults.\n", filename);
        return;
    }

    char line[MAX_LINE];
    char current_section[64] = "";
    char current_label[64] = "";
    int current_design = -1;
    
    // 首先设置默认值
    init_default_config(config);

    while (fgets(line, sizeof(line), file)) {
        // 跳过空行和注释
        if (line[0] == ';' || line[0] == '#' || line[0] == '\n') 
            continue;
            
        // 移除行尾的换行符
        line[strcspn(line, "\r\n")] = 0;
        
        // 检查是否是节段标记 [Section]
        if (line[0] == '[' && line[strlen(line)-1] == ']') {
            strncpy(current_section, line + 1, strlen(line) - 2);
            current_section[strlen(line) - 2] = '\0';
            continue;
        }
        
        // 解析键值对
        char* key = strtok(line, "=");
        char* value = strtok(NULL, "\n");
        
        if (!key || !value) continue;
        
        key = trim(key);
        value = trim(value);
        
        // 根据当前节段处理配置
        if (str_equals(current_section, "General")) {
            if (str_equals(key, "k"))
                config->k = atoi(value);
            else if (str_equals(key, "threads"))
                config->threads = atoi(value);
            else if (str_equals(key, "min_overlap"))
                config->min_overlap = atoi(value);
            else if (str_equals(key, "max_mismatch_ratio"))
                config->max_mismatch_ratio = atof(value);
        }
        else if (str_equals(current_section, "Input")) {
            if (str_equals(key, "fq1"))
                strncpy(config->fq1, value, sizeof(config->fq1) - 1);
            else if (str_equals(key, "fq2"))
                strncpy(config->fq2, value, sizeof(config->fq2) - 1);
        }
        else if (str_equals(current_section, "Output")) {
            if (str_equals(key, "output_file"))
                strncpy(config->output_file, value, sizeof(config->output_file) - 1);
            else if (str_equals(key, "output_format"))
                config->output_format = atoi(value);
        }
        else if (str_equals(current_section, "Advanced")) {
            if (str_equals(key, "batch_size"))
                config->batch_size = atoi(value);
            else if (str_equals(key, "quality_threshold"))
                config->quality_threshold = atoi(value);
        }
        else if (str_equals(current_section, "Design")) {
            if (str_equals(key, "seq")) {
                // Convert sequence to uppercase before adding design
                char* uppercase_seq = strdup(value);
                for(int i = 0; uppercase_seq[i]; i++) {
                    uppercase_seq[i] = toupper(uppercase_seq[i]);
                }
                add_design(&config->designs, uppercase_seq, current_label);
                free(uppercase_seq);
                current_design = config->designs.count - 1;
            }
            else if (str_equals(key, "label")) {
                strncpy(current_label, value, 63);
                current_label[63] = '\0';
            }
        }
        else if (str_equals(current_section, "ZONE") && current_design >= 0) {
            static char zone_name[64] = "";
            static int zone_start = -1;
            static int zone_end = -1;
            
            if (str_equals(key, "name"))
                strncpy(zone_name, value, 63);
            else if (str_equals(key, "start"))
                zone_start = atoi(value);
            else if (str_equals(key, "end")) {
                zone_end = atoi(value);
                
                if (zone_name[0] != '\0' && zone_start != -1 && zone_end != -1) {
                    add_zone(&config->designs.designs[current_design].zones, 
                            zone_name, zone_start, zone_end,
                            config->designs.designs[current_design].seq,
                            config->k);  // 传入序列和k值
                    
                    zone_name[0] = '\0';
                    zone_start = -1;
                    zone_end = -1;
                }
            }
        }
    }

    fclose(file);
}

// Function to display configurations
void print_config(config_t* config) {
    printf("\n=== Configuration Settings ===\n");
    
    printf("\n[General]\n");
    printf("k-mer size: %d\n", config->k);
    printf("Threads: %d\n", config->threads);
    printf("Minimum overlap: %d\n", config->min_overlap);
    printf("Maximum mismatch ratio: %.2f\n", config->max_mismatch_ratio);
    
    printf("\n[Input]\n");
    printf("Read1 file: %s\n", config->fq1);
    printf("Read2 file: %s\n", config->fq2[0] ? config->fq2 : "None");
    printf("Mode: %s\n", config->fq2[0] ? "Paired-end" : "Single-end");
    
    printf("\n[Output]\n");
    printf("Output file: %s\n", config->output_file);
    printf("Output format: %d\n", config->output_format);
    
    printf("\n[Advanced]\n");
    printf("Batch size: %d\n", config->batch_size);
    printf("Quality threshold: %d\n", config->quality_threshold);
    
    printf("\n[Designs]\n");
    printf("Number of designs: %d\n", config->designs.count);
    
    for (int d = 0; d < config->designs.count; d++) {
        printf("\nDesign %d:\n", d + 1);
        design_t* design = &config->designs.designs[d];
        
        printf("Label: %s\n", design->label);  // 打印label
        if (design->seq) {
            size_t seq_len = strlen(design->seq);
            printf("Sequence length: %zu\n", seq_len);
            printf("Sequence:\n");
            for (size_t i = 0; i < seq_len; i++) {
                printf("%c", design->seq[i]);
                if ((i + 1) % 60 == 0) printf("\n");
            }
            if (seq_len % 60 != 0) printf("\n");
        }
        
        printf("Number of zones: %d\n", design->zones.count);
        for (int i = 0; i < design->zones.count; i++) {
            printf("Zone %d: %s (positions %d-%d)\n", 
                   i + 1,
                   design->zones.zones[i].name,
                   design->zones.zones[i].start,
                   design->zones.zones[i].end);
            printf("  Left k-mer: %s\n", design->zones.zones[i].left_kmer);
            printf("  Right k-mer: %s\n", design->zones.zones[i].right_kmer);
        }
    }
    
    printf("\n==============================\n\n");
}



// Initialize kmer set
void init_kmer_set(kmer_set_t *kmer_set) {
    kmer_set->set = kh_init(str);
}

// Compute kmer set from a sequence
void compute_kmer_set(const char *seq, int k, kmer_set_t *kmer_set) {
    int seq_len = strlen(seq);
    if (seq_len < k) return;  // Add length check
    
    for (int i = 0; i <= seq_len - k; i++) {
        char *kmer = (char *)malloc(k + 1);
        strncpy(kmer, seq + i, k);
        kmer[k] = '\0';
        
        // Convert k-mer to uppercase for consistent comparison
        for (int j = 0; j < k; j++) {
            kmer[j] = toupper(kmer[j]);
        }
        
        int ret;
        khiter_t iter = kh_put(str, kmer_set->set, kmer, &ret);
        if (!ret) {
            // The kmer was already in the set
            free(kmer);
        }
    }
}

// Check if kmer is in set
int kmer_in_set(kmer_set_t *kmer_set, const char *kmer) {
    // Create temporary uppercase version of kmer for comparison
    char *upper_kmer = strdup(kmer);
    for (int i = 0; upper_kmer[i]; i++) {
        upper_kmer[i] = toupper(upper_kmer[i]);
    }
    
    khiter_t iter = kh_get(str, kmer_set->set, upper_kmer);
    int result = (iter != kh_end(kmer_set->set));
    
    free(upper_kmer);
    return result;
}

// Free kmer set
void free_kmer_set(kmer_set_t *kmer_set) {
    const char *key;
    khiter_t k;
    for (k = kh_begin(kmer_set->set); k != kh_end(kmer_set->set); ++k) {
        if (kh_exist(kmer_set->set, k)) {
            key = kh_key(kmer_set->set, k);
            free((char *)key);
        }
    }
    kh_destroy(str, kmer_set->set);
}

// Function to compute Jaccard index
float compute_jaccard_index(kmer_set_t *set1, kmer_set_t *set2) {
    int intersection = 0;
    int union_count = 0;
    
    // Debug output
    printf("Set sizes - Set1: %d, Set2: %d\n", kh_size(set1->set), kh_size(set2->set));

    // Count intersection and first part of union
    for (khiter_t k = kh_begin(set1->set); k != kh_end(set1->set); ++k) {
        if (!kh_exist(set1->set, k)) continue;
        
        const char *kmer = kh_key(set1->set, k);
        if (kmer_in_set(set2, kmer)) {
            intersection++;
        }
        union_count++;
    }

    // Add remaining elements from set2 to union
    for (khiter_t k = kh_begin(set2->set); k != kh_end(set2->set); ++k) {
        if (!kh_exist(set2->set, k)) continue;
        
        const char *kmer = kh_key(set2->set, k);
        if (!kmer_in_set(set1, kmer)) {
            union_count++;
        }
    }

    // Debug output
    printf("Intersection: %d, Union: %d\n", intersection, union_count);

    if (union_count == 0) return 0.0;
    return (float)intersection / union_count;
}

// Function to concatenate two reads with reverse complement of seq2 and 100 A's spacer
char* concatenate_reads(const char* seq1, const char* seq2) {
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    const int spacer_len = 100;  // Length of poly-A spacer
    
    // Allocate memory for concatenated sequence (seq1 + spacer + rev_comp_seq2)
    char* result = (char*)malloc(len1 + spacer_len + len2 + 1);
    
    // Copy seq1
    strcpy(result, seq1);
    
    // Add poly-A spacer
    for(int i = 0; i < spacer_len; i++) {
        result[len1 + i] = 'A';
    }
    
    // Create reverse complement of seq2
    char* rev_comp_seq2 = (char*)malloc(len2 + 1);
    reverse_complement(seq2, rev_comp_seq2, len2);
    
    // Add reverse complemented seq2
    strcpy(result + len1 + spacer_len, rev_comp_seq2);
    
    // Free temporary memory
    free(rev_comp_seq2);
    
    result[len1 + spacer_len + len2] = '\0';
    return result;
}

// Function to find kmer position in a sequence
int find_kmer_position(const char* seq, int seq_len, const char* kmer, int k, int start_pos) {
    for (int i = start_pos; i <= seq_len - k; i++) {
        if (strncmp(seq + i, kmer, k) == 0) {
            return i;
        }
    }
    return -1;
}

// Modified worker function
void process_reads(void *_data, long i, int tid) {
    process_batch_t *data = (process_batch_t*)_data;
    read_pair_t *pair = &data->pairs[i];

    char *assembled_seq = NULL;
    int assembled_length = 0;

    if (data->is_pe) {
        // Process paired-end reads
        char merged_seq[MAX_READ_LEN * 2];
        char merged_qual[MAX_READ_LEN * 2];
        assembled_length = pe_assemble(pair->r1.seq, pair->r1.qual,
                                       pair->r2.seq, pair->r2.qual,
                                       data->min_overlap, data->max_mismatch_ratio,
                                       merged_seq, merged_qual);
        if (assembled_length > 0) {
            assembled_seq = strdup(merged_seq);
        } else {
            // If unable to assemble, concatenate reads
            assembled_seq = concatenate_reads(pair->r1.seq, pair->r2.seq);
            assembled_length = strlen(assembled_seq);
        }
    } else {
        // Process single-end reads
        assembled_seq = strdup(pair->r1.seq);
        assembled_length = strlen(assembled_seq);
    }

    // Determine the best matching design using Jaccard Index
    int best_design_index = -1;
    float best_jaccard = 0.0;
    kmer_set_t read_kmers;
    
    // Print assembled sequence
    printf("Assembled sequence: %s\n", assembled_seq);
    printf("Assembled length: %d\n", (int)strlen(assembled_seq));
    
    init_kmer_set(&read_kmers);
    compute_kmer_set(assembled_seq, data->k, &read_kmers);

    printf("k-mer size: %d\n", data->k);

    printf("\nComparing with designs:\n");

    // printf("\nRead k-mers:\n");
    // for (khiter_t k = kh_begin(read_kmers.set); k != kh_end(read_kmers.set); ++k) {
    //     if (kh_exist(read_kmers.set, k)) {
    //         printf("  %s\n", kh_key(read_kmers.set, k));
    //     }
    // }
    printf("Number of designs: %d\n", data->designs->count);

    for (int d = 0; d < data->designs->count; d++) {
        float jaccard = compute_jaccard_index(&read_kmers, &data->designs->designs[d].kmer_set);

        // printf("\nDesign k-mers:\n");
        // for (khiter_t k = kh_begin(data->designs->designs[d].kmer_set.set); k != kh_end(data->designs->designs[d].kmer_set.set); ++k) {
        //     if (kh_exist(data->designs->designs[d].kmer_set.set, k)) {
        //         printf("  %s\n", kh_key(data->designs->designs[d].kmer_set.set, k));
        //     }
        // }
        
        printf("Design %d (%s):\n", d, data->designs->designs[d].label);
        printf("  - Sequence: %s\n", data->designs->designs[d].seq);
        printf("  - Jaccard Index: %.4f\n", jaccard);
        
        if (jaccard > best_jaccard) {
            best_jaccard = jaccard;
            best_design_index = d;
            printf("  -> New best match!\n");
        }
    }
    
    if (best_design_index >= 0) {
        printf("\nBest matching design: %d (%s)\n", 
               best_design_index, 
               data->designs->designs[best_design_index].label);
        printf("Final Jaccard Index: %.4f\n", best_jaccard);
    } else {
        printf("\nNo matching design found\n");
    }

    if (best_design_index >= 0) {
        // Process the zones of the best matching design
        design_t *best_design = &data->designs->designs[best_design_index];
        int is_reverse = 0;
        int zone_found = 0;

        // First try on the original assembled sequence
        zone_found = process_zones(assembled_seq, assembled_length, data->k, best_design);

        if (!zone_found) {
            // Try reverse complement
            char *rev_comp_seq = (char *)malloc((assembled_length + 1) * sizeof(char));
            reverse_complement(assembled_seq, rev_comp_seq, assembled_length);
            is_reverse = 1;
            zone_found = process_zones(rev_comp_seq, assembled_length, data->k, best_design);
            free(rev_comp_seq);
        }

        // Output the results
        output_results(pair, best_design, is_reverse);
    } else {
        // No matching design found
        // Handle as needed (e.g., output a message or skip)
    }

    // Clean up
    free(assembled_seq);
    free_kmer_set(&read_kmers);
}

// Function to reverse complement a DNA sequence
void reverse_complement(const char* seq, char* rev_comp_seq, int length) {
    for (int i = 0; i < length; i++) {
        char base = seq[length - i - 1];
        switch(base) {
            case 'A': rev_comp_seq[i] = 'T'; break;
            case 'C': rev_comp_seq[i] = 'G'; break;
            case 'G': rev_comp_seq[i] = 'C'; break;
            case 'T': rev_comp_seq[i] = 'A'; break;
            case 'N': rev_comp_seq[i] = 'N'; break;
            default:  rev_comp_seq[i] = 'N'; break;
        }
    }
    rev_comp_seq[length] = '\0';
}

// Function to reverse a quality string
void reverse_quality(const char* qual, char* rev_qual, int length) {
    for (int i = 0; i < length; i++) {
        rev_qual[i] = qual[length - i - 1];
    }
    rev_qual[length] = '\0';
}

// Function to merge two reads (replacement for assemble_pairs)
int pe_assemble(const char* read1_seq, const char* read1_qual,
                const char* read2_seq, const char* read2_qual,
                int min_overlap, float max_mismatch_fraction,
                char* merged_seq, char* merged_qual) {
    // Initialize reads
    int read1_length = strlen(read1_seq);
    int read2_length = strlen(read2_seq);

    // Reverse complement read2
    char* rev_comp_seq = (char*)malloc((read2_length + 1) * sizeof(char));
    char* rev_qual = (char*)malloc((read2_length + 1) * sizeof(char));
    reverse_complement(read2_seq, rev_comp_seq, read2_length);
    reverse_quality(read2_qual, rev_qual, read2_length);

    // Prepare for alignment
    int max_overlap = (read1_length < read2_length) ? read1_length : read2_length;
    int best_overlap = 0;
    float best_mismatch_fraction = 1.0;
    int best_position = -1;

    // Try all possible overlaps
    for (int overlap = max_overlap; overlap >= min_overlap; overlap--) {
        int mismatches = 0;
        int matches = 0;

        // Compare the overlapping region
        for (int i = 0; i < overlap; i++) {
            char base1 = read1_seq[read1_length - overlap + i];
            char base2 = rev_comp_seq[i];
            if (base1 != base2 && base1 != 'N' && base2 != 'N') {
                mismatches++;
            } else {
                matches++;
            }
        }

        float mismatch_fraction = (float)mismatches / (matches + mismatches);
        if (mismatch_fraction <= max_mismatch_fraction) {
            if (mismatch_fraction < best_mismatch_fraction ||
                (mismatch_fraction == best_mismatch_fraction && overlap > best_overlap)) {
                best_overlap = overlap;
                best_mismatch_fraction = mismatch_fraction;
                best_position = read1_length - overlap;
            }
        }
    }

    if (best_position == -1) {
        // No valid overlap found
        free(rev_comp_seq);
        free(rev_qual);
        return 0; // Indicate failure
    } else {
        // Merge the reads
        int prefix_length = best_position;
        int suffix_length = read2_length - best_overlap;

        // Copy prefix from read1
        strncpy(merged_seq, read1_seq, prefix_length);
        // Merge overlapping region
        for (int i = 0; i < best_overlap; i++) {
            char base1 = read1_seq[prefix_length + i];
            char base2 = rev_comp_seq[i];
            char qual1 = read1_qual[prefix_length + i];
            char qual2 = rev_qual[i];

            // Choose base with higher quality
            if (qual1 >= qual2) {
                merged_seq[prefix_length + i] = base1;
                merged_qual[prefix_length + i] = qual1;
            } else {
                merged_seq[prefix_length + i] = base2;
                merged_qual[prefix_length + i] = qual2;
            }
        }
        // Copy suffix from reversed read2
        strncpy(merged_seq + prefix_length + best_overlap, rev_comp_seq + best_overlap, suffix_length);
        strncpy(merged_qual + prefix_length + best_overlap, rev_qual + best_overlap, suffix_length);

        // Null-terminate the strings
        int merged_length = prefix_length + best_overlap + suffix_length;
        merged_seq[merged_length] = '\0';
        merged_qual[merged_length] = '\0';

        free(rev_comp_seq);
        free(rev_qual);
        return merged_length; // Return the length of merged read
    }
}

// Function to initialize read structure
void init_read(read_t *r) {
    r->name = NULL;
    r->seq = NULL;
    r->qual = NULL;
}

// Function to free read structure
void free_read(read_t *r) {
    free(r->name);
    free(r->seq);
    free(r->qual);
}

// Function to process zones of a design
int process_zones(const char *assembled_seq, int assembled_length, int k, design_t *design) {
    int zone_found = 0;

    for (int i = 0; i < design->zones.count; i++) {
        zone_t *zone = &design->zones.zones[i];
        int left_pos = find_kmer_position(assembled_seq, assembled_length, zone->left_kmer, k, 0);
        int right_pos = -1;
        char left_found = '-';
        char right_found = '-';

        if (left_pos >= 0) {
            left_found = '+';
            int zone_seq_start = left_pos + k;
            int zone_seq_end = zone_seq_start + (zone->end - zone->start);
            if (zone_seq_end <= assembled_length) {
                // Extract zone sequence
                char *zone_seq = (char *)malloc((zone_seq_end - zone_seq_start + 1) * sizeof(char));
                strncpy(zone_seq, assembled_seq + zone_seq_start, zone_seq_end - zone_seq_start);
                zone_seq[zone_seq_end - zone_seq_start] = '\0';

                // Check right k-mer
                right_pos = find_kmer_position(assembled_seq, assembled_length, zone->right_kmer, k, zone_seq_end);
                if (right_pos == zone_seq_end) {
                    right_found = '+';
                }

                // Output result
                printf("%s\t%s\t%s\t%c\t%c\t", design->label, zone->name, zone_seq, left_found, right_found);
                free(zone_seq);
                zone_found = 1;
            }
        } else {
            // Try finding right k-mer
            right_pos = find_kmer_position(assembled_seq, assembled_length, zone->right_kmer, k, 0);
            if (right_pos >= 0) {
                right_found = '+';
                int zone_seq_end = right_pos;
                int zone_seq_start = zone_seq_end - (zone->end - zone->start);
                if (zone_seq_start >= 0) {
                    // Extract zone sequence
                    char *zone_seq = (char *)malloc((zone_seq_end - zone_seq_start + 1) * sizeof(char));
                    strncpy(zone_seq, assembled_seq + zone_seq_start, zone_seq_end - zone_seq_start);
                    zone_seq[zone_seq_end - zone_seq_start] = '\0';

                    // Output result
                    printf("%s\t%s\t%s\t%c\t%c\t", design->label, zone->name, zone_seq, left_found, right_found);
                    free(zone_seq);
                    zone_found = 1;
                }
            }
        }
    }

    printf("\n");
    return zone_found;
}

// Function to output results
void output_results(read_pair_t *pair, design_t *design, int is_reverse) {
    // Output the results as needed
    // This function can be customized based on your output requirements
    // Example placeholder output:
    printf("Read: %s\n", pair->r1.name);
    printf("Design: %s\n", design->label);
    printf("Is reverse: %d\n", is_reverse);
}

int main(int argc, char *argv[]) {
    int c;
    char *fq1 = NULL, *fq2 = NULL;
    char config_file[1024] = "config.ini";  // Default config file
    config_t config;

    // Parse command line arguments
    while ((c = getopt(argc, argv, "1:2:k:t:c:")) >= 0) {
        switch (c) {
            case '1': fq1 = optarg; break;
            case '2': fq2 = optarg; break;
            case 'k': g_k = atoi(optarg); break;
            case 't': n_threads = atoi(optarg); break;
            case 'c': strncpy(config_file, optarg, sizeof(config_file) - 1); break;
        }
    }

    // Load config file
    load_config(config_file, &config);

    // Command line arguments override config file
    if (!fq1 && config.fq1[0] != '\0') fq1 = config.fq1;
    if (!fq2 && config.fq2[0] != '\0') fq2 = config.fq2;
    if (g_k == 31) g_k = config.k;  // Only override if not set in command line
    if (n_threads == DEFAULT_THREADS) n_threads = config.threads;

    if (!fq1) {
        fprintf(stderr, "Usage: %s [-c config.ini] [-1 <read1.fq>] [-2 <read2.fq>] [-k kmer_size] [-t threads]\n", argv[0]);
        return 1;
    }

    // Print final configuration
    print_config(&config);

    // Open FASTQ files
    gzFile fp1 = gzopen(fq1, "r");
    if (!fp1) {
        fprintf(stderr, "Error opening file %s\n", fq1);
        return 1;
    }
    gzFile fp2 = fq2 ? gzopen(fq2, "r") : NULL;
    if (fq2 && !fp2) {
        fprintf(stderr, "Error opening file %s\n", fq2);
        return 1;
    }
    kseq_t *ks1 = kseq_init(fp1);
    kseq_t *ks2 = fq2 ? kseq_init(fp2) : NULL;

    // Process reads in batches
    read_pair_t *pairs = malloc(DEFAULT_BLOCK_SIZE * sizeof(read_pair_t));
    process_batch_t batch = {
        .pairs = pairs,
        .n_pairs = 0,
        .k = g_k,
        .is_pe = fq2 != NULL,
        .min_overlap = config.min_overlap,
        .max_mismatch_ratio = config.max_mismatch_ratio,
        .designs = &config.designs  // Pass the designs
    };

    while (1) {
        batch.n_pairs = 0;

        // Read a batch of sequences
        while (batch.n_pairs < DEFAULT_BLOCK_SIZE) {
            if (kseq_read(ks1) < 0) break;

            init_read(&pairs[batch.n_pairs].r1);
            pairs[batch.n_pairs].r1.name = strdup(ks1->name.s);
            pairs[batch.n_pairs].r1.seq = strdup(ks1->seq.s);
            pairs[batch.n_pairs].r1.qual = strdup(ks1->qual.s);

            if (fq2) {
                if (kseq_read(ks2) < 0) break;
                init_read(&pairs[batch.n_pairs].r2);
                pairs[batch.n_pairs].r2.name = strdup(ks2->name.s);
                pairs[batch.n_pairs].r2.seq = strdup(ks2->seq.s);
                pairs[batch.n_pairs].r2.qual = strdup(ks2->qual.s);
            }

            batch.n_pairs++;
        }

        if (batch.n_pairs == 0) break;

        // Process the batch
        kt_for(n_threads, process_reads, &batch, batch.n_pairs);

        // Clean up batch
        for (int i = 0; i < batch.n_pairs; i++) {
            free_read(&pairs[i].r1);
            if (fq2) {
                free_read(&pairs[i].r2);
            }
        }
    }

    // Clean up
    free(pairs);
    kseq_destroy(ks1);
    if (ks2) kseq_destroy(ks2);
    gzclose(fp1);
    if (fp2) gzclose(fp2);

    // Clean up designs
    free_designs_array(&config.designs);

    return 0;
}

char* extract_kmer(const char* seq, int position, int k, int direction) {
    char* kmer = (char*)malloc((k + 1) * sizeof(char));
    if (direction == -1) {
        // Extract k bases to the left
        if (position - k < 0) {
            free(kmer);
            return NULL;
        }
        strncpy(kmer, seq + position - k, k);
    } else {
        // Extract k bases to the right
        int seq_len = strlen(seq);
        if (position + k > seq_len) {
            free(kmer);
            return NULL;
        }
        strncpy(kmer, seq + position, k);
    }
    kmer[k] = '\0';
    return kmer;
}
