/* Deterministic sanitizer fallback; this is not coverage-guided fuzzing. */
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
int LLVMFuzzerTestOneInput(const uint8_t *, size_t);
int main(int argc,char **argv) {
    uint32_t state=argc>2?(uint32_t)strtoul(argv[2],NULL,10):20261002u;
    size_t runs=argc>1?(size_t)strtoul(argv[1],NULL,10):10000;
    uint8_t input[128];
    for(size_t i=0;i<runs;i++) {
        for(size_t j=0;j<sizeof(input);j++) { state=state*1664525u+1013904223u;input[j]=(uint8_t)(state>>24); }
        LLVMFuzzerTestOneInput(input,4+i%125);
    }
    return 0;
}
