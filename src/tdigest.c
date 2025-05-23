#include "tdigest.h"
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef TD_MALLOC_INCLUDE
#define TD_MALLOC_INCLUDE "td_malloc.h"
#endif

#include TD_MALLOC_INCLUDE

#define __td_max(x, y) (((x) > (y)) ? (x) : (y))
#define __td_min(x, y) (((x) < (y)) ? (x) : (y))

#define TD_CAPACITY_MULTIPLIER 6
#define TD_CAPACITY_FIXEDOFFSET 10
#define CDF_MEDIAN 0.5
#define PICO 1e-12
#define HALF 0.5

static inline double weighted_average_sorted(double x1, double w1, double x2, double w2) {
    const double x = (x1 * w1 + x2 * w2) / (w1 + w2);
    return __td_max(x1, __td_min(x, x2));
}

static inline bool tdigest_long_long_add_safe(long long a, long long b) {
    if (b < 0) {
        return (a >= __LONG_LONG_MAX__ - b);
    }
    return (a <= __LONG_LONG_MAX__ - b);
}

static inline double weighted_average(double x1, double w1, double x2, double w2) {
    if (x1 <= x2) {
        return weighted_average_sorted(x1, w1, x2, w2);
    }
    return weighted_average_sorted(x2, w2, x1, w1);
}

static inline void swap(double *arr, unsigned int i, unsigned int j) {
    const double temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

static inline void swap_l(long long *arr, unsigned int i, unsigned int j) {
    const long long temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

static unsigned int partition(double *means, long long *weights, unsigned int start,
                              unsigned int end, unsigned int pivot_idx) {
    const double pivotMean = means[pivot_idx];
    swap(means, pivot_idx, end);
    swap_l(weights, pivot_idx, end);

    unsigned int i = start - 1;

    for (unsigned int j = start; j < end; j++) {
        // If current element is smaller than the pivot
        if (means[j] < pivotMean) {
            // increment index of smaller element
            i++;
            swap(means, i, j);
            swap_l(weights, i, j);
        }
    }
    swap(means, i + 1, end);
    swap_l(weights, i + 1, end);
    return i + 1;
}

/**
 * Standard quick sort except that sorting rearranges parallel arrays
 *
 * @param means  Values to sort on
 * @param weights The auxillary values to sort.
 * @param start  The beginning of the values to sort
 * @param end    The value after the last value to sort
 */
static void td_qsort(double *means, long long *weights, unsigned int start, unsigned int end) {
    if (start < end) {
        // two elements can be directly compared
        if ((end - start) == 1) {
            if (means[start] > means[end]) {
                swap(means, start, end);
                swap_l(weights, start, end);
            }
            return;
        }
        // generating a random number as a pivot was very expensive vs the array size
        // const unsigned int pivot_idx = start + rand()%(end - start + 1);
        const unsigned int pivot_idx = (end + start) / 2; // central pivot
        const unsigned int new_pivot_idx = partition(means, weights, start, end, pivot_idx);
        if (new_pivot_idx > start) {
            td_qsort(means, weights, start, new_pivot_idx - 1);
        }
        td_qsort(means, weights, new_pivot_idx + 1, end);
    }
}

static inline size_t cap_from_compression(double compression) {
    if ((size_t)compression >
        ((SIZE_MAX / sizeof(double) / TD_CAPACITY_MULTIPLIER) - TD_CAPACITY_FIXEDOFFSET)) {
        return 0;
    }
    return (TD_CAPACITY_MULTIPLIER * (size_t)(compression)) + TD_CAPACITY_FIXEDOFFSET;
}

static inline bool should_td_compress(td_histogram_t *h) {
    return ((h->merged_nodes + h->unmerged_nodes) >= (h->cap - 1));
}

static inline unsigned long next_node(td_histogram_t *h) {
    return h->merged_nodes + h->unmerged_nodes;
}

static inline int check_overflow(const double v) {
    // double-precision overflow detected on h->unmerged_weight
    if (isinf(v)) {
        return EDOM;
    }
    return 0;
}

static inline int check_td_overflow(const double new_unmerged_weight,
                                    const double new_total_weight) {
    // double-precision overflow detected on h->unmerged_weight
    if (isinf(new_unmerged_weight)) {
        return EDOM;
    }
    if (isinf(new_total_weight)) {
        return EDOM;
    }
    const double denom = 2 * MM_PI * new_total_weight * log(new_total_weight);
    if (isinf(denom)) {
        return EDOM;
    }

    return 0;
}

unsigned long td_centroid_count(td_histogram_t *h) { return next_node(h); }

void td_reset(td_histogram_t *h) {
    if (!h) {
        return;
    }
    h->min = __DBL_MAX__;
    h->max = -h->min;
    h->merged_nodes = 0;
    h->merged_weight = 0;
    h->unmerged_nodes = 0;
    h->unmerged_weight = 0;
    h->total_compressions = 0;
}

int td_init(double compression, td_histogram_t **result) {

    const size_t capacity = cap_from_compression(compression);
    if (capacity < 1) {
        return 1;
    }
    td_histogram_t *histogram = NULL;
    histogram = (td_histogram_t *)__td_malloc(sizeof(td_histogram_t));
    if (!histogram) {
        return 1;
    }
    // set them to null to avoid uninitialized branch in td_free
    histogram->nodes_mean = NULL;
    histogram->nodes_weight = NULL;
    histogram->cap = capacity;
    histogram->compression = (double)compression;
    td_reset(histogram);
    histogram->nodes_mean = (double *)__td_calloc(capacity, sizeof(double));
    if (!histogram->nodes_mean) {
        td_free(histogram);
        return 1;
    }
    histogram->nodes_weight = (long long *)__td_calloc(capacity, sizeof(long long));
    if (!histogram->nodes_weight) {
        td_free(histogram);
        return 1;
    }
    *result = histogram;

    return 0;
}

td_histogram_t *td_new(double compression) {
    td_histogram_t *mdigest = NULL;
    td_init(compression, &mdigest);
    return mdigest;
}

void td_free(td_histogram_t *histogram) {
    if (histogram->nodes_mean) {
        __td_free((void *)(histogram->nodes_mean));
    }
    if (histogram->nodes_weight) {
        __td_free((void *)(histogram->nodes_weight));
    }
    __td_free((void *)(histogram));
}

int td_merge(td_histogram_t *into, td_histogram_t *from) {
    if (td_compress(into) != 0) {
        return EDOM;
    }
    if (td_compress(from) != 0) {
        return EDOM;
    }
    const unsigned long pos = from->merged_nodes + from->unmerged_nodes;
    for (unsigned long i = 0; i < pos; i++) {
        const double mean = from->nodes_mean[i];
        const long long weight = from->nodes_weight[i];
        if (td_add(into, mean, weight) != 0) {
            return EDOM;
        }
    }
    return 0;
}

long long td_size(td_histogram_t *h) { return h->merged_weight + h->unmerged_weight; }

double td_cdf(td_histogram_t *h, double val) {
    td_compress(h);
    // no data to examine
    if (h->merged_nodes == 0) {
        return (double)NAN;
    }
    // below lower bound
    if (val < h->min) {
        return 0;
    }
    // above upper bound
    if (val > h->max) {
        return 1;
    }
    if (h->merged_nodes == 1) {
        // exactly one centroid, should have max==min
        const double width = h->max - h->min;
        if (val - h->min <= width) {
            // min and max are too close together to do any viable interpolation
            return CDF_MEDIAN;
        }
        // interpolate if somehow we have width > 0 and max != min
        return (val - h->min) / width;
    }
    const unsigned long n = h->merged_nodes;
    // check for the left tail
    const double left_centroid_mean = h->nodes_mean[0];
    const double left_centroid_weight = (double)h->nodes_weight[0];
    const double merged_weight_d = (double)h->merged_weight;
    if (val < left_centroid_mean) {
        // note that this is different than h->nodes_mean[0] > min
        // ... this guarantees we divide by non-zero number and interpolation works
        const double width = left_centroid_mean - h->min;
        if (width > 0) {
            // must be a sample exactly at min
            if (fabs(val - h->min) < PICO) {
                return CDF_MEDIAN / merged_weight_d;
            }
            return (1 + (val - h->min) / width * (left_centroid_weight / 2 - 1)) / merged_weight_d;
        }
        // this should be redundant with the check val < h->min
        return 0;
    }
    // and the right tail
    const double right_centroid_mean = h->nodes_mean[n - 1];
    const double right_centroid_weight = (double)h->nodes_weight[n - 1];
    if (val > right_centroid_mean) {
        const double width = h->max - right_centroid_mean;
        if (width > 0) {
            if (fabs(val - h->max) < PICO) {
                return 1 - CDF_MEDIAN / merged_weight_d;
            }
            // there has to be a single sample exactly at max
            const double dq =
                (1 + (h->max - val) / width * (right_centroid_weight / 2 - 1)) / merged_weight_d;
            return 1 - dq;
        }
        return 1;
    }
    // we know that there are at least two centroids and mean[0] < x < mean[n-1]
    // that means that there are either one or more consecutive centroids all at exactly x
    // or there are consecutive centroids, c0 < x < c1
    double weightSoFar = 0;
    for (unsigned long it = 0; it < n - 1; it++) {
        // weightSoFar does not include weight[it] yet
        if (fabs(h->nodes_mean[it] - val) < PICO) {
            // we have one or more centroids == x, treat them as one
            // dw will accumulate the weight of all of the centroids at x
            double dw = 0;
            while (it < n && fabs(h->nodes_mean[it] - val) < PICO) {
                dw += (double)h->nodes_weight[it];
                it++;
            }
            return (weightSoFar + dw / 2) / (double)h->merged_weight;
        }
        if (h->nodes_mean[it] <= val && val < h->nodes_mean[it + 1]) {
            const double node_weight = (double)h->nodes_weight[it];
            const double node_weight_next = (double)h->nodes_weight[it + 1];
            const double node_mean = h->nodes_mean[it];
            const double node_mean_next = h->nodes_mean[it + 1];
            // landed between centroids ... check for floating point madness
            if (node_mean_next - node_mean > 0) {
                // note how we handle singleton centroids here
                // the point is that for singleton centroids, we know that their entire
                // weight is exactly at the centroid and thus shouldn't be involved in
                // interpolation
                double leftExcludedW = 0;
                double rightExcludedW = 0;
                if (fabs(node_weight - 1) < PICO) {
                    if (fabs(node_weight_next - 1) < PICO) {
                        // two singletons means no interpolation
                        // left singleton is in, right is out
                        return (weightSoFar + 1) / merged_weight_d;
                    }
                    leftExcludedW = HALF;
                } else if (fabs(node_weight_next - 1) < PICO) {
                    rightExcludedW = HALF;
                }
                double dw = (node_weight + node_weight_next) / 2;

                // adjust endpoints for any singleton
                double dwNoSingleton = dw - leftExcludedW - rightExcludedW;

                double base = weightSoFar + node_weight / 2 + leftExcludedW;
                return (base + dwNoSingleton * (val - node_mean) / (node_mean_next - node_mean)) /
                       merged_weight_d;
            } // this is simply caution against floating point madness
              // it is conceivable that the centroids will be different
              // but too near to allow safe interpolation
            double dw = (node_weight + node_weight_next) / 2;
            return (weightSoFar + dw) / merged_weight_d;
        }
        weightSoFar += (double)h->nodes_weight[it];
    }
    return 1 - CDF_MEDIAN / merged_weight_d;
}

static double td_internal_iterate_centroids_to_index(const td_histogram_t *h, const double index,
                                                     const double left_centroid_weight,
                                                     const unsigned long total_centroids,
                                                     double *weightSoFar, unsigned long *node_pos) {
    if (left_centroid_weight > 1 && index < left_centroid_weight / 2) {
        // there is a single sample at min so we interpolate with less weight
        return h->min + (index - 1) / (left_centroid_weight / 2 - 1) * (h->nodes_mean[0] - h->min);
    }

    // usually the last centroid will have unit weight so this test will make it moot
    if (index > (double)h->merged_weight - 1) {
        return h->max;
    }

    // if the right-most centroid has more than one sample, we still know
    // that one sample occurred at max so we can do some interpolation
    const double right_centroid_weight = (double)h->nodes_weight[total_centroids - 1];
    const double right_centroid_mean = h->nodes_mean[total_centroids - 1];
    if (right_centroid_weight > 1 &&
        (double)h->merged_weight - index <= right_centroid_weight / 2) {
        return h->max - ((double)h->merged_weight - index - 1) / (right_centroid_weight / 2 - 1) *
                            (h->max - right_centroid_mean);
    }

    for (; *node_pos < total_centroids - 1; (*node_pos)++) {
        const unsigned long i = *node_pos;
        const double node_weight = (double)h->nodes_weight[i];
        const double node_weight_next = (double)h->nodes_weight[i + 1];
        const double node_mean = h->nodes_mean[i];
        const double node_mean_next = h->nodes_mean[i + 1];
        const double dw = (node_weight + node_weight_next) / 2;
        if (*weightSoFar + dw > index) {
            // centroids i and i+1 bracket our current point
            // check for unit weight
            double leftUnit = 0;
            if (fabs(node_weight - 1.0) < PICO) {
                if (index - *weightSoFar < HALF) {
                    // within the singleton's sphere
                    return node_mean;
                }
                leftUnit = HALF;
            }
            double rightUnit = 0;
            if (fabs(node_weight_next - 1.0) < PICO) {
                if (*weightSoFar + dw - index <= HALF) {
                    // no interpolation needed near singleton
                    return node_mean_next;
                }
                rightUnit = HALF;
            }
            const double z1 = index - *weightSoFar - leftUnit;
            const double z2 = *weightSoFar + dw - index - rightUnit;
            return weighted_average(node_mean, z2, node_mean_next, z1);
        }
        *weightSoFar += dw;
    }

    // weightSoFar = totalWeight - weight[total_centroids-1]/2 (very nearly)
    // so we interpolate out to max value ever seen
    const double z1 = index - (double)h->merged_weight - right_centroid_weight / 2.0;
    const double z2 = right_centroid_weight / 2 - z1;
    return weighted_average(right_centroid_mean, z1, h->max, z2);
}

double td_quantile(td_histogram_t *h, double q) {
    td_compress(h);
    // q should be in [0,1]
    if (q < 0.0 || q > 1.0 || h->merged_nodes + h->unmerged_nodes == 0) {
        return (double)NAN;
    }
    // with one data point, all quantiles lead to Rome
    if (h->merged_nodes + h->unmerged_nodes == 1) {
        return h->nodes_mean[0];
    }

    // if values were stored in a sorted array, index would be the offset we are interested in
    const double index = q * (double)h->merged_weight;

    // beyond the boundaries, we return min or max
    // usually, the first centroid will have unit weight so this will make it moot
    if (index < 1) {
        return h->min;
    }

    // we know that there are at least two centroids now
    const unsigned long n = h->merged_nodes;

    // if the left centroid has more than one sample, we still know
    // that one sample occurred at min so we can do some interpolation
    const double left_centroid_weight = (double)h->nodes_weight[0];

    // in between extremes we interpolate between centroids
    double weightSoFar = left_centroid_weight / 2;
    unsigned long i = 0;
    return td_internal_iterate_centroids_to_index(h, index, left_centroid_weight, n, &weightSoFar,
                                                  &i);
}

int td_quantiles(td_histogram_t *h, const double *quantiles, double *values, size_t length) {
    td_compress(h);

    if (NULL == quantiles || NULL == values) {
        return EINVAL;
    }

    const unsigned long n = h->merged_nodes;
    if (n == 0) {
        for (size_t i = 0; i < length; i++) {
            values[i] = (double)NAN;
        }
        return 0;
    }
    if (n == 1) {
        for (size_t i = 0; i < length; i++) {
            const double requested_quantile = quantiles[i];

            // q should be in [0,1]
            if (requested_quantile < 0.0 || requested_quantile > 1.0) {
                values[i] = (double)NAN;
            } else {
                // with one data point, all quantiles lead to Rome
                values[i] = h->nodes_mean[0];
            }
        }
        return 0;
    }

    // we know that there are at least two centroids now
    // if the left centroid has more than one sample, we still know
    // that one sample occurred at min so we can do some interpolation
    const double left_centroid_weight = (double)h->nodes_weight[0];

    // in between extremes we interpolate between centroids
    double weightSoFar = left_centroid_weight / 2;
    unsigned long node_pos = 0;

    // to avoid allocations we use the values array for intermediate computation
    // i.e. to store the expected cumulative count at each percentile
    for (size_t qpos = 0; qpos < length; qpos++) {
        const double index = quantiles[qpos] * (double)h->merged_weight;
        values[qpos] = td_internal_iterate_centroids_to_index(h, index, left_centroid_weight, n,
                                                              &weightSoFar, &node_pos);
    }
    return 0;
}

static double td_internal_trimmed_mean(const td_histogram_t *h, const double leftmost_weight,
                                       const double rightmost_weight) {
    double count_done = 0;
    double trimmed_sum = 0;
    double trimmed_count = 0;
    for (unsigned long i = 0; i < h->merged_nodes; i++) {

        const double n_weight = (double)h->nodes_weight[i];
        // Assume the whole centroid falls into the range
        double count_add = n_weight;

        // If we haven't reached the low threshold yet, skip appropriate part of the centroid.
        count_add -= __td_min(__td_max(0, leftmost_weight - count_done), count_add);

        // If we have reached the upper threshold, ignore the overflowing part of the centroid.

        count_add = __td_min(__td_max(0, rightmost_weight - count_done), count_add);

        // consider the whole centroid processed
        count_done += n_weight;

        // increment the sum / count
        trimmed_sum += h->nodes_mean[i] * count_add;
        trimmed_count += count_add;

        // break once we cross the high threshold
        if (count_done >= rightmost_weight) {
            break;
        }
    }

    return trimmed_sum / trimmed_count;
}

double td_trimmed_mean_symmetric(td_histogram_t *h, double proportion_to_cut) {
    td_compress(h);
    // proportion_to_cut should be in [0,1]
    if (h->merged_nodes == 0 || proportion_to_cut < 0.0 || proportion_to_cut > 1.0) {
        return (double)NAN;
    }
    // with one data point, all values lead to Rome
    if (h->merged_nodes == 1) {
        return h->nodes_mean[0];
    }

    /* translate the percentiles to counts */
    const double leftmost_weight = floor((double)h->merged_weight * proportion_to_cut);
    const double rightmost_weight = ceil((double)h->merged_weight * (1.0 - proportion_to_cut));

    return td_internal_trimmed_mean(h, leftmost_weight, rightmost_weight);
}

double td_trimmed_mean(td_histogram_t *h, double leftmost_cut, double rightmost_cut) {
    td_compress(h);
    // leftmost_cut and rightmost_cut should be in [0,1]
    if (h->merged_nodes == 0 || leftmost_cut < 0.0 || leftmost_cut > 1.0 || rightmost_cut < 0.0 ||
        rightmost_cut > 1.0) {
        return (double)NAN;
    }
    // with one data point, all values lead to Rome
    if (h->merged_nodes == 1) {
        return h->nodes_mean[0];
    }

    /* translate the percentiles to counts */
    const double leftmost_weight = floor((double)h->merged_weight * leftmost_cut);
    const double rightmost_weight = ceil((double)h->merged_weight * rightmost_cut);

    return td_internal_trimmed_mean(h, leftmost_weight, rightmost_weight);
}

int td_add(td_histogram_t *h, double val, long long weight) {
    if (should_td_compress(h)) {
        const int overflow_res = td_compress(h);
        if (overflow_res != 0) {
            return overflow_res;
        }
    }
    const unsigned long pos = next_node(h);
    if (pos >= h->cap) {
        return EDOM;
    }
    if (tdigest_long_long_add_safe(h->unmerged_weight, weight) == false) {
        return EDOM;
    }
    const long long new_unmerged_weight = h->unmerged_weight + weight;
    if (tdigest_long_long_add_safe(new_unmerged_weight, h->merged_weight) == false) {
        return EDOM;
    }
    const long long new_total_weight = new_unmerged_weight + h->merged_weight;
    // double-precision overflow detected
    const int overflow_res =
        check_td_overflow((double)new_unmerged_weight, (double)new_total_weight);
    if (overflow_res != 0) {
        return overflow_res;
    }

    if (val < h->min) {
        h->min = val;
    }
    if (val > h->max) {
        h->max = val;
    }
    h->nodes_mean[pos] = val;
    h->nodes_weight[pos] = weight;
    h->unmerged_nodes++;
    h->unmerged_weight = new_unmerged_weight;
    return 0;
}

int td_compress(td_histogram_t *h) {
    if (h->unmerged_nodes == 0) {
        return 0;
    }
    unsigned long N = h->merged_nodes + h->unmerged_nodes;
    td_qsort(h->nodes_mean, h->nodes_weight, 0, (unsigned int)N - 1);
    const double total_weight = (double)h->merged_weight + (double)h->unmerged_weight;
    // double-precision overflow detected
    const int overflow_res = check_td_overflow((double)h->unmerged_weight, (double)total_weight);
    if (overflow_res != 0) {
        return overflow_res;
    }
    if (total_weight <= 1) {
        return 0;
    }
    const double denom = 2 * MM_PI * total_weight * log(total_weight);
    if (check_overflow(denom) != 0) {
        return EDOM;
    }

    // Compute the normalizer given compression and number of points.
    const double normalizer = h->compression / denom;
    if (check_overflow(normalizer) != 0) {
        return EDOM;
    }
    unsigned long cur = 0;
    double weight_so_far = 0;

    for (unsigned long i = 1; i < N; i++) {
        const double proposed_weight = (double)h->nodes_weight[cur] + (double)h->nodes_weight[i];
        const double z = proposed_weight * normalizer;
        // quantile up to cur
        const double q0 = weight_so_far / total_weight;
        // quantile up to cur + i
        const double q2 = (weight_so_far + proposed_weight) / total_weight;
        // Convert  a quantile to the k-scale
        const bool should_add = (z <= (q0 * (1 - q0))) && (z <= (q2 * (1 - q2)));
        // next point will fit
        // so merge into existing centroid
        if (should_add) {
            h->nodes_weight[cur] += h->nodes_weight[i];
            const double delta = h->nodes_mean[i] - h->nodes_mean[cur];
            const double weighted_delta =
                (delta * (double)h->nodes_weight[i]) / (double)h->nodes_weight[cur];
            h->nodes_mean[cur] += weighted_delta;
        } else {
            weight_so_far += (double)h->nodes_weight[cur];
            cur++;
            h->nodes_weight[cur] = h->nodes_weight[i];
            h->nodes_mean[cur] = h->nodes_mean[i];
        }
        if (cur != i) {
            h->nodes_weight[i] = 0;
            h->nodes_mean[i] = 0.0;
        }
    }
    h->merged_nodes = cur + 1;
    h->merged_weight = (long long)total_weight;
    h->unmerged_nodes = 0;
    h->unmerged_weight = 0;
    h->total_compressions++;
    return 0;
}

double td_min(td_histogram_t *h) { return h->min; }

double td_max(td_histogram_t *h) { return h->max; }

double td_compression(td_histogram_t *h) { return h->compression; }

const long long *td_centroids_weight(td_histogram_t *h) { return h->nodes_weight; }

const double *td_centroids_mean(td_histogram_t *h) { return h->nodes_mean; }

long long td_centroids_weight_at(td_histogram_t *h, unsigned long pos) {
    return h->nodes_weight[pos];
}

double td_centroids_mean_at(td_histogram_t *h, unsigned long pos) {
    if (pos > h->merged_nodes) {
        return (double)NAN;
    }
    return h->nodes_mean[pos];
}
