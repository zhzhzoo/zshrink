#include "zshrink.h"
#include "zshrink_utilities.h"
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void
display_help()
{
    printf("zshrink -- compress an alnum file.\n"
           "\tusage zshrink x <input> <output> : extract a file\n"
           "\t      zshrink c <input> <output> [-e] : compress a file\n"
           "\t                                        -e to enable encryption\n");
}

int
main(int argc, char *argv[])
{
    zshrink *f;
    char *pwd = 0;
    int ret;
    if (argc < 4) {
        display_help();
        return 0;
    }

    if (strcmp(argv[1], "x") == 0) {
        DEBUGPRINT("extracting\n");
        f = zshrink_open(argv[2]);
        if (!f) {
            printf("Input file doesn't exist or isn't a valid zsr archive...\n");
            exit(0);
        }
        if (zshrink_encrypted(f)) {
            pwd = "1995415";//getpass("It's encrypted.\nPassword:");
        }
        ret = zshrink_decompress(f, argv[3], pwd);
        if (ret != 0) {
            printf("Some error occured..\n");
            return ret;
        }
        zshrink_close(f);
    }

    else if(strcmp(argv[1], "c") == 0) {
        DEBUGPRINT("compressing\n");
        pwd = 0;
        if (argc >= 5 && strcmp(argv[4], "-e") == 0) {
            pwd = "1995415";//getpass("Password:");
        }
        f = zshrink_compress(argv[2], argv[3], pwd);
        if (!f) {
            printf("Some error occured..\n");
            return -1;
        }
        zshrink_close(f);
    }

    else {
        display_help();
        return 0;
    }

    return 0;
}
