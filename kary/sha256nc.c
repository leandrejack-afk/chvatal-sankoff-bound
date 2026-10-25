/* sha256nc FILE...: SHA-256 of large files, read with F_NOCACHE so that they do not fill the
 * page cache (CommonCrypto; output in the format of shasum -a 256). */
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <CommonCrypto/CommonDigest.h>
int main(int argc, char **argv) {
    const size_t B = (size_t)1 << 26;
    unsigned char *buf = malloc(B);
    for (int i = 1; i < argc; i++) {
        int fd = open(argv[i], O_RDONLY);
        if (fd < 0) { perror(argv[i]); return 1; }
        fcntl(fd, F_NOCACHE, 1);
        CC_SHA256_CTX c;
        CC_SHA256_Init(&c);
        ssize_t n;
        while ((n = read(fd, buf, B)) > 0) CC_SHA256_Update(&c, buf, (CC_LONG)n);
        if (n < 0) { perror("read"); return 1; }
        close(fd);
        unsigned char d[CC_SHA256_DIGEST_LENGTH];
        CC_SHA256_Final(d, &c);
        for (int j = 0; j < CC_SHA256_DIGEST_LENGTH; j++) printf("%02x", d[j]);
        printf("  %s\n", argv[i]);
    }
    return 0;
}
