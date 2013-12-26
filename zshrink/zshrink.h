#ifndef HAS_ZSHRINK_H
#include <stdio.h>

/* zsr 压缩文件格式
 * 字节
 * 1-2      0x89 0x64 文件标识
 * 3-4      0x00 不加密 0x01 加密
 * 5-8      原文件CRC32
 * 9-16     压缩前长度
 * 17-24    压缩后长度
 * 25-+oo   数据
 */

typedef struct {
    char                *file_name;
    FILE                *fp;
    unsigned int        encrypted;
    unsigned int        crc32;
    unsigned long long  length_before;
    unsigned long long  length_after;
} zshrink;

zshrink *zshrink_open(char *name);
int zshrink_decompress(zshrink *zs, char *dest_name, char *key);
void zshrink_close(zshrink *zs);

zshrink *zshrink_compress(char *src_name, char *dest_name, char *key);
unsigned int zshrink_encrypted(zshrink *zs);

#endif
