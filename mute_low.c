#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "wavefile.h"


static const int kernel_len = 767; //255; /* Must be odd */
static double high_pass_kernel[1023];
static double low_pass_kernel[1023];
static double bandpass_kernel[1023];

static const char prefix[] = "ml_";
#define FILTER_PERCENTAGE 18
#define CUTOFF_HIGH (1000.0)  /* Frequency in Hz */
#define CUTOFF_LOW  (  50.0)  /* Frequency in Hz */
#define PI (3.141592653589793238462643383)

#define max(a,b)  ((a)>(b) ? (a) : (b))

#define PATH_SEPERATOR '/'

void calc_filter(double cutoff_high, double cutoff_low) {
    int i;
    double total = 0;

    for(i = 0; i < kernel_len; i++) {
        int index = i -(kernel_len-1)/2;
        double angle;
        angle = index *2.0*PI*cutoff_high;
        if(index == 0)
            low_pass_kernel[i] = 1.0;
        else
            low_pass_kernel[i] = sin(angle)/angle;
    }
    for(i = 0; i < kernel_len; i++) {
        int index = i-(kernel_len-1)/2;
        low_pass_kernel[i] *= ( 0.42+0.5*cos(2.0*PI*index/(kernel_len-1))+0.08*cos(4.0*PI*index/(kernel_len-1)));
    }
    total = 0.0;
    for(i = 0; i < kernel_len; i++) {
        total += low_pass_kernel[i];
    }
    for(i = 0; i < kernel_len; i++) {
        low_pass_kernel[i] /= total;
    }

    for(i = 0; i < kernel_len; i++) {
        int index = i -(kernel_len-1)/2;
        double angle;
        angle = index *2.0*PI*cutoff_low;
        if(index == 0)
            high_pass_kernel[i] = 1.0;
        else
            high_pass_kernel[i] = sin(angle)/angle;
    }
    for(i = 0; i < kernel_len; i++) {
        int index = i-(kernel_len-1)/2;
        high_pass_kernel[i] *= ( 0.42+0.5*cos(2.0*PI*index/(kernel_len-1))+0.08*cos(4.0*PI*index/(kernel_len-1)));
    }
    total = 0.0;
    for(i = 0; i < kernel_len; i++) {
        total += high_pass_kernel[i];
    }
    for(i = 0; i < kernel_len; i++) {
        high_pass_kernel[i] /= total;
    }

    for(i = 0; i < kernel_len; i++) {
        int index = i-(kernel_len-1)/2;
        bandpass_kernel[i] = low_pass_kernel[i] - high_pass_kernel[i];
    }
}

struct wave *filter_bandpass(struct wave *s) {
    struct wave *n;
    int i, new_len;
    double *totals;
    totals = malloc(sizeof(double) * s->channel_count); 
    if(totals == NULL) {
      return NULL;
    }

    /* Genearte a new file name by appending the prefix */

    /* Work out how long the new file will be, and allocate it */
    new_len = s->sample_count-(kernel_len-1);
    n = wavefile_new(s->sample_rate,new_len, s->channel_count);
    if(n == NULL) {
       free(totals);
       return NULL;
    }

    /* Now filter and decimate the input samples */
    for(i = 0; i < new_len; i++) {
        int c,n_channels;
        n_channels = s->channel_count;
        for(c = 0; c < n_channels; c++) {
           totals[c] = 0.0;
        }

        /* Filter the input */
        for(c = 0; c < n_channels; c++) {
          int j;
        
          for(j = 0; j < kernel_len; j++) {
              totals[c] += (double)s->channel_data[c][i+j]  * bandpass_kernel[j];
          }
          n->channel_data[c][i]  = (int)(totals[c]);
        }
        if(i%(s->sample_rate*5)==0)
            printf("%i:%02i seconds processed\n",i/(s->sample_rate)/60, i/(s->sample_rate)%60);
    }
    free(totals);
    return n;
}

struct wave *filter_rumble(struct wave *s) {
    struct wave *n;
    int i, new_len;
    double *totals;
    totals = malloc(sizeof(double) * s->channel_count);
    if(totals == NULL) {
      return NULL;
    }
    /* Genearte a new file name by appending the prefix */

    /* Work out how long the new file will be, and allocate it */
    new_len = s->sample_count-(kernel_len-1);
    n = wavefile_new(s->sample_rate,new_len,s->channel_count);
    if(n == NULL) {
       free(totals);
       return NULL;
    }

    /* Now filter and decimate the input samples */
    for(i = 0; i < new_len; i++) {
        int c, n_channels = s->channel_count;
        for(c = 0; c < n_channels; c++) {
           totals[c] = 0.0;
        }
#if 0
        int j;

        /* Filter the input */
        for(j = 0; j < kernel_len; j++) {
            total_l += (double)s->left_samples[i+j]  * high_pass_kernel[j];
            total_r += (double)s->right_samples[i+j] * high_pass_kernel[j];
        }
        n->left_samples[i]  = (int)(s->left_samples[i+kernel_len/2]  - total_l);
        n->right_samples[i] = (int)(s->right_samples[i+kernel_len/2] - total_r);

        if(i%(s->sample_rate*5)==0)
            printf("%i:%02i seconds processed\n",i/(s->sample_rate)/60, i/(s->sample_rate)%60);
#else
        for(c = 0; c < n_channels; c++) {
          n->channel_data[c][i]  = (int)(s->channel_data[c][i+kernel_len/2]  - totals[c]);
        }
#endif
    }
    free(totals);
    return n;
}

void power_per_ms(struct wave *a,struct wave *b) {
    int i, percent;
    int s_per_ms = b->sample_rate/10;
    int blocks = (b->sample_count-1)/s_per_ms;
    int *levels;
    int *filtered_levels;
    long long muted = 0,unmuted = 0;

    levels = malloc(blocks*sizeof(int));
    if(levels == NULL) {
      printf("Out of memory in power_per_ms()\n");
      return;
    }

    filtered_levels = malloc(blocks*sizeof(int));
    if(filtered_levels == NULL) {
      free(levels);
      printf("Out of memory in power_per_ms()\n");
      return;
    }

    for(i=0; i < blocks; i++)
    {
        int j;
        unsigned long long total_a = 0;
        unsigned long long total_b = 0; 
        int c, n_channels = a->channel_count;
        for(c = 0; c < n_channels; c++) { 
           for(j = 0; j < s_per_ms; j++)
           {
            
              int s_a, s_b;
              s_a = a->channel_data[c][i*s_per_ms+j];
              total_a += s_a*s_a;
              s_b = b->channel_data[c][i*s_per_ms+j];
              total_b += s_b*s_b;
           }
        }
        if(total_a < 40000*kernel_len)
           total_a = 40000*kernel_len;
        levels[i] = total_a > 99 ? total_b/(total_a/100) : 0;

        if(levels[i]>FILTER_PERCENTAGE)
          printf("%i:%02i.%1i:  %10lli  %10lli   %10i\n", i/600,i/10%60,i%10, total_a, total_b, levels[i]);
    }
#if 0
    for(i = 0; i < blocks; i++) {
        filtered_levels[i] = 0;

        if(i > 3) filtered_levels[i] = max(filtered_levels[i], levels[i-4]);
        if(i > 2) filtered_levels[i] = max(filtered_levels[i], levels[i-3]);
        if(i > 1) filtered_levels[i] = max(filtered_levels[i], levels[i-2]);
        if(i > 0) filtered_levels[i] = max(filtered_levels[i], levels[i-1]);
                          filtered_levels[i] = max(filtered_levels[i], levels[i]);
        if(i < blocks-1) filtered_levels[i] = max(filtered_levels[i], levels[i+1]);
        if(i < blocks-2) filtered_levels[i] = max(filtered_levels[i], levels[i+2]);
        if(i < blocks-3) filtered_levels[i] = max(filtered_levels[i], levels[i+3]);
        if(i < blocks-4) filtered_levels[i] = max(filtered_levels[i], levels[i+4]);
        if(i < blocks-5) filtered_levels[i] = max(filtered_levels[i], levels[i+5]);
        if(i < blocks-6) filtered_levels[i] = max(filtered_levels[i], levels[i+6]);
    }
    for(i = 0; i < a->sample_count; i++) {
       int block = i / s_per_ms;
       if(block < 0)        block = 0;
       if(block > blocks-1) block = blocks-1;
       if(filtered_levels[block]>FILTER_PERCENTAGE) {
           a->left_samples[i]  = 0;
           a->right_samples[i] = 0;
        muted++;
       } else {
        unmuted++;
       }
    }
#else
    for(i = 0; i < blocks; i++) {
        int j;
        filtered_levels[i] = 0;
        for(j=-5; j < 10;j++)
        {
            if(i+j >= 0 && i+j< blocks) if(levels[i+j]>FILTER_PERCENTAGE) filtered_levels[i]++;
        }
    }
    for(i = 0; i < a->sample_count; i++) {
       int block = (i- kernel_len/2) / s_per_ms;
       if(block < 0)        block = 0;
       if(block > blocks-1) block = blocks-1;
       if(filtered_levels[block]>0) {
           int c, n_channels = a->channel_count;
           for(c = 0; c < n_channels; c++)
             a->channel_data[c][i]  = 0;
        muted++;
       } else {
        unmuted++;
       }
    }

#endif
    percent = muted*100/(unmuted+muted);
    printf("%i%% muted\n",percent);

    free(filtered_levels);
    free(levels);
}

int main(int argc, char *argv[]) {
    int rtn = 0;
    while(argc > 1) {
        struct wave *s;
        s = wavefile_read(argv[1]);
        if(s != NULL) {
            struct wave *n_bp,*n_hp;
            calc_filter(CUTOFF_HIGH/(float)s->sample_rate, CUTOFF_LOW/(float)s->sample_rate);

            n_bp = filter_bandpass(s);
            n_hp = filter_rumble(s);
            /* Mute based on power per block in low freq */
            if(n_bp != NULL && n_hp != NULL)
                power_per_ms(n_hp,n_bp);

            if(n_bp != NULL && n_hp != NULL) {
                char *new_name;
                int name_len = strlen(argv[1]);
                int prefix_len = strlen(prefix);
                new_name = malloc(name_len+prefix_len+1);
                if(new_name != NULL) {
                    int i,j;

                    /* Copy the file name over */
                    strcpy(new_name,argv[1]);

                    /* Shuffle the filename forward until either we are at the 
                     * start of the name, or we hit the directory seperator */
                    i = name_len;
                    new_name[i+prefix_len] = new_name[i];
                    while(i > 0 && new_name[i-1] != PATH_SEPERATOR) {
                      i--; 
                      new_name[i+prefix_len] = new_name[i];
                    }

                    /* Merge in the prefix */
                    for(j = 0; j < prefix_len; j++)
                       new_name[i+j] = prefix[j];

                    if(wavefile_write(n_hp, new_name)) {
                       printf("Wrote output as '%s'\n",new_name);
                    } else {
                       printf("Unable to write output as '%s'\n",new_name);
                       rtn = 1;
                    }
                } else {
                  printf("Out of memory\n");
                  rtn = 1;
                }
            }
            if(n_bp != NULL)    wavefile_destroy(n_bp);
            if(n_hp != NULL)    wavefile_destroy(n_hp);
            wavefile_destroy(s);
        } else {
            printf("Unable to read file %s\n",argv[1]);
            rtn = 1;
        }
        argc--;
        argv++;
    }
    return rtn;
}
