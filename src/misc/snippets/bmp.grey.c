/* write a grayscale bitmap */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned char  BYTE;  /* one byte (0-255) */
typedef unsigned short WORD;  /* two bytes */
typedef unsigned int   DWORD; /* four bytes */

/* Bitmap format cf.
 * http://en.wikipedia.org/wiki/BMP_file_format
 * http://upload.wikimedia.org/wikipedia/commons/c/c4/BMPfileFormat.png
 * http://www.fileformat.info/format/bmp/egff.htm */

typedef struct {
/* BYTE type[2]; */       /* = "BM", it is manipulated separately to make
                             the size of this structure a multiple of 4 */
  DWORD sizeFile;         /* = total file size == offset + bitmap-size */
  DWORD reserved;         /* == 0 */
  DWORD offset;           /* offset from start of file == sizeof(BitmapHeader) with "BM"
                             + sizeof(DIBHeader) + size of palette */
} BitmapHeader;

typedef struct {
  DWORD sizeStruct;       /* sizeof(DIBHeader) */
  DWORD width, height;
  WORD  planes;           /* 1 */
  WORD  bitCount;         /* bits of each pixel, 8 for 256-color, 24 for RGB true color */
  DWORD compression;      /* 0 */
  DWORD sizeImage;        /* (width+?)(till multiple of 4) * height in bytes */
  DWORD xPixelsPerMeter;  /* resolution in mspaint, 2952 */
  DWORD yPixelsPerMeter;  /* resolution in mspaint, 2952 */
  DWORD colorUsed;        /* 256 for 256-color, 0 for true color */
  DWORD colorImportant;   /* 256 for 256-color, 0 for true color */
#if 0
  DWORD maskRed;          /* 0x00ff0000 */
  DWORD maskGreen;        /* 0x0000ff00 */
  DWORD maskBlue;         /* 0x000000ff */
  DWORD maskAlpha;        /* 0xff000000 */
  DWORD colorSpaceType;   /* 1 */
  DWORD colorSpaceEnd[9]; /* 0, 0, 0xfc1e5486, 0, 0, 0xfc666669, 0, 0, 0xff28f5c4 */
  DWORD gammaRed;         /* 0 */
  DWORD gammaGreen;       /* 0 */
  DWORD gammaBlue;        /* 0 */
  DWORD Intent;           /* 4 */
  DWORD ICCData;          /* 0 */
  DWORD ICCSize;          /* 0 */
  DWORD Reserve;          /* 0 */
#endif
} DIBHeader;


typedef struct {
  size_t w;
  size_t ww;  /* ceil4(width) */
  size_t h;
  size_t size;
  unsigned char *data;
} bmp_t;

#define RGB(r, g, b) ((((r)&0xFF)<<16) + (((g)&0xFF)<<8) + ((b)&0xFF))
#define CEIL4(x) ((((x)+3)/4)*4)

static void bmp_close(bmp_t *b) {
  if (b->data) free(b->data);
  free(b);
}

/* create a new bitmap file */
static bmp_t *bmp_open(int width, int height)
{
  bmp_t *b;

  if ((b = calloc(1, sizeof(*b))) == NULL) exit(1);

  b->w = (size_t) width;
  b->h = (size_t) height;
  b->ww = CEIL4(b->w);
  b->size = b->ww * b->h;

  if ((b->data = calloc(b->size, 1)) == NULL) exit(1);
  return b;
}

static int bmp_save(const bmp_t *b, const char *fn)
{
  FILE *fp;
  BitmapHeader bh = {0};
  DIBHeader dh = {0};
  DWORD pal[256];
  int i;

  bh.offset = sizeof(BitmapHeader) + sizeof(DIBHeader) + sizeof(pal) + 2;
  bh.sizeFile = bh.offset + b->size;

  dh.sizeStruct = sizeof(DIBHeader);
  dh.width = b->w;
  dh.height = b->h;
  dh.planes = 1;
  dh.sizeImage = b->size;
  dh.bitCount = 8;
  dh.colorUsed = dh.colorImportant = 256;
  printf("file %u, struct %u\n", bh.sizeFile, dh.sizeStruct);

  for (i = 0; i < 256; i++) {
    pal[255 - i] = RGB(i, i, i);
  }

  if ((fp = fopen(fn, "wb")) == NULL) {
    fprintf(stderr, "cannot save to BMP file %s\n", fn);
    return 1;
  }
  fputc('B', fp);
  fputc('M', fp); /* write type */
  /* fill BMP file header */
  fwrite(&bh, sizeof bh, 1, fp);
  fwrite(&dh, sizeof dh, 1, fp);
  fwrite(pal, sizeof pal, 1, fp);
  if (fwrite(b->data, 1, b->size, fp) != b->size) {
    fclose(fp);
    return 1;
  }
  fclose(fp);
  return 0;
}

__inline static void bmp_dot(bmp_t *b, DWORD x, DWORD y, DWORD color)
{
  if (x < b->w && y < b->h)
    *(b->data + b->ww * y + x) = (BYTE) color;
}


static void writebmp(const char *fn)
{
  int i;
  bmp_t *bmp;

  bmp = bmp_open(16, 16);
  for (i = 0; i < 16; i++) bmp_dot(bmp, 3, i, 125);
  bmp_save(bmp, fn);
  bmp_close(bmp);
}


int main(void)
{
  writebmp("out.bmp");
  return 0;
}

