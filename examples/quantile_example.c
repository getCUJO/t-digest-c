#include <stdio.h>

#include "tdigest.h"

#define STREAM_SIZE 1000000

static inline double randMToN(double M, double N)
{
	//NOLINTNEXTLINE(cert-msc30-c, cert-msc50-cpp)
	return M + (rand() / (RAND_MAX / (N - M)));
}

int main(void)
{
	td_histogram_t *mdigest =
		td_new(500); //NOLINT(readability-magic-numbers)
	printf("compression is %f capacity is %lu\n",
	       mdigest->compression,
	       mdigest->cap);
	double seeds[STREAM_SIZE];
	for (int i = 0; i < STREAM_SIZE; ++i) {
		seeds[i] = randMToN(0, 10);
	}

	for (int i = 0; i < STREAM_SIZE; ++i) {
		td_add(mdigest, seeds[i], 1);
	}
	td_compress(mdigest);
	for (int i = 0; i < 10; ++i) {
		const double v = seeds[i];
		printf("value %f is at percentile %f\n", v, td_cdf(mdigest, v));
	}
	printf("\n");
	for (int i = 0; i <= 100; i += 10) { //NOLINT(readability-magic-numbers)
		printf("%d percentile has value %f\n",
		       i,
		       td_quantile(mdigest, i / 100.0));
	}
}
