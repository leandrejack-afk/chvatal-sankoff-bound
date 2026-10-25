/* prints kbig's inputs count and weight types for (k, h) without building: make_inputs only */
#define main kbig_main
#include "../../kbig.c"
#undef main
int main(int argc, char **argv) {
    K = atoi(argv[1]); H = atoi(argv[2]);
    M1 = (1ULL << H) - 1; M2 = (1ULL << (2 * H)) - 1;
    make_inputs();
    fprintf(stderr, "k=%d h=%d inputs=%d weight types=%d\n", K, H, NI, NW);
    for (int i = 0; i < NW; i++)
        printf("%d %d %d %d %d %d\n", i, wtypes[i].size, wtypes[i].ta, wtypes[i].tb, wtypes[i].na, wtypes[i].nb);
    return 0;
}
