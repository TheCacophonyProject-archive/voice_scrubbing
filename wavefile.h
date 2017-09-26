/******************************************************************************
* An API for reading or writing WAV files
******************************************************************************/

/* Details of the memory structure used to hold the WAV file data */
struct wave {
    int total_size;
    int channel_count;
    int sample_rate;
    int data_rate;
    int block_align;
    int bits_per_sample;
    int sample_count;
    int *left_samples;
    int *right_samples;
};

/* Definitions of the three functions */
struct wave *wavefile_new(int sample_rate, int sample_count);
struct wave *wavefile_read(char *filename);
int wavefile_write(struct wave *w, char *filename);
void wavefile_destroy(struct wave *w);
