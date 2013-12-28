#include "zshrink.h"
#include "zshrink_utilities.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "zomb_datastream.h"
#include "zomb_arithmetic_coder.h"
#include "zomb_fileio.h"
#include "zomb_encrypt.h"

zshrink
*zshrink_open(char *name) {
    zshrink *zsr;
    FILE *fp;
    unsigned char type[4];

    DEBUGPRINT("Opening existing file %s\n", name);
    fp = fopen(name, "rb");
    if (fp == 0) {
        DEBUGPRINT("File open failed\n");
        return 0;
    }

    zsr = malloc(sizeof(zshrink));

    fread(type, 1, 4, fp);
    if (type[0] != 0x89 || type[1] != 0x64 || type[2] != 0) {
        DEBUGPRINT("File type not supported\n");
        free(zsr);
        fclose(fp);
        return 0;
    }

    if (type[3] == 0x00) {
        zsr->encrypted = 0;
    }
    else if (type[3] == 0x01) {
        zsr->encrypted = 1;
    }
    else {
        DEBUGPRINT("File type not supported\n");
        free(zsr);
        fclose(fp);
        return 0;
    }

    fread(&zsr->crc32, 4, 1, fp);
    fread(&zsr->length_before, 8, 1, fp);
    fread(&zsr->length_after, 8, 1, fp);

    DEBUGPRINT("Successful\n"
               "\tencrypted : %d crc32 : %x"
               "\tlength_before : %llu length_after : %llu\n",
               zsr->encrypted, zsr->crc32, zsr->length_before, zsr->length_after);
    zsr->fp = fp;

    zsr->file_name = malloc((strlen(name) + 1) * sizeof(char));
    memcpy(zsr->file_name, name, strlen(name));

    return zsr;
}

int /* 0 success, -1 fail */
zshrink_decompress(zshrink *zs, char *dest_name, char *key)
{
    FILE *fout;
    zomb_data_processor_t no_encrypt[3];
    zomb_data_processor_t encrypt[4];
    void *arg[4], **ret;

    no_encrypt[0] = zomb_filein;
    no_encrypt[1] = zomb_decoder;
    no_encrypt[2] = zomb_fileout;
    encrypt[0] = zomb_filein;
    encrypt[1] = zomb_encrypt;
    encrypt[2] = zomb_decoder;
    encrypt[3] = zomb_fileout;
    DEBUGPRINT("Opening output %s\n", dest_name);
    fout = fopen(dest_name, "wb");
    if (fout == 0) {
        DEBUGPRINT("... failed%s\n", dest_name);
        return -1;
    }

    if (zs->encrypted) {
        if (key == 0 || key[0] == 0) {
            return -1;
        }

        fseek(zs->fp, 24, SEEK_SET);
        arg[0] = zs->fp;
        arg[1] = malloc(sizeof(char) * strlen(key));
        memcpy(arg[1], key, strlen(key));
        arg[2] = &(zs->length_before);
        arg[3] = fout;

        DEBUGPRINT("Extracting encrypted\n");
        zomb_datastream_run(encrypt, 4, arg, &ret);

        fclose(fout);
        free(arg[1]);
        free(ret);
        return 0;
    }
    else {
        fseek(zs->fp, 24, SEEK_SET);
        arg[0] = zs->fp;
        arg[1] = &(zs->length_before);
        arg[2] = fout;

        DEBUGPRINT("Extracting non-encrypted\n");
        zomb_datastream_run(no_encrypt, 3, arg, &ret);

        fclose(fout);
        free(ret);
        return 0;
    }
}

void
zshrink_close(zshrink *zs) {
    DEBUGPRINT("Closing zshrink\n");
    fclose(zs->fp);
    free(zs->file_name);
    free(zs);
}

zshrink*
zshrink_compress(char *src_name, char *dest_name, char *key) {
    FILE *fin, *fout;
    zomb_data_processor_t no_encrypt[3];
    zomb_data_processor_t encrypt[4];
    zshrink *zs;
    void *arg[4], **ret;
    unsigned char zero = 0, encrypted;
    unsigned char type[2] = {0x89, 0x64};
    int i;

    no_encrypt[0] = zomb_filein;
    no_encrypt[1] = zomb_encoder;
    no_encrypt[2] = zomb_fileout;
    encrypt[0] = zomb_filein;
    encrypt[1] = zomb_encoder;
    encrypt[2] = zomb_encrypt;
    encrypt[3] = zomb_fileout;
    DEBUGPRINT("Opening input %s\n", src_name);
    fin = fopen(src_name, "rb");
    if (fin == 0) {
        DEBUGPRINT("... failed\n");
        return 0;
    }

    DEBUGPRINT("Opening output %s\n", dest_name);
    fout = fopen(dest_name, "wb");
    if (fout == 0) {
        DEBUGPRINT("... failed\n");
        return 0;
    }
    fwrite(&zero, 1, 24, fout);

    zs = malloc(sizeof(zshrink));
    if (key) {
        arg[0] = fin;
        arg[1] = 0;
        arg[2] = malloc(sizeof(char) * strlen(key));
        memcpy(arg[2], key, strlen(key));
        arg[3] = fout;

        DEBUGPRINT("Compressing encrypted\n");
        zomb_datastream_run(encrypt, 4, arg, &ret);

        zs->fp = fout;
        zs->length_before = *((long long*)ret[0]);
        zs->length_after = *((long long*)ret[1]);

        free(arg[2]);
        free(ret[0]);
        free(ret[1]);
        free(ret[2]);
        free(ret[3]);

        encrypted = 1;
        zs->encrypted = 1;
        zs->crc32 = 0;
        zs->file_name = malloc(sizeof(char) * strlen(dest_name));
        fclose(fin);
    }
    else {
        arg[0] = fin;
        arg[1] = 0;
        arg[2] = fout;

        DEBUGPRINT("Compressing unencrypted\n");
        zomb_datastream_run(no_encrypt, 3, arg, &ret);

        zs->fp = fout;
        zs->length_before = *((long long*)ret[0]);
        zs->length_after = *((long long*)ret[1]);

        free(ret[0]);
        free(ret[1]);
        free(ret[2]);

        encrypted = 0;
        zs->encrypted = 0;
        zs->crc32 = 0;
        zs->file_name = malloc(sizeof(char) * strlen(dest_name));
        fclose(fin);
    }

    DEBUGPRINT("Writing metas\n");
    DEBUGPRINT("\tencrypted : %d crc32 : %x"
               "\tlength_before : %llu length_after : %llu\n",
               zs->encrypted, zs->crc32, zs->length_before, zs->length_after);
    rewind(zs->fp);
    fwrite(type, 1, 2, zs->fp);
    fwrite(&zero, 1, 1, zs->fp);
    fwrite(&encrypted, 1, 1, zs->fp);
    for (i = 0; i < 4; i++) {
        fwrite(&zero, 1, 1, zs->fp);
    }
    fwrite(&zs->length_before, 8, 1, zs->fp);
    fwrite(&zs->length_after, 8, 1, zs->fp);

    return zs;
}

unsigned int /* returns a boolean indicating whether it's encrypted */
zshrink_encrypted(zshrink *zs)
{
    return zs->encrypted;
}
