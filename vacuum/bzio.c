/*

*/
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define WORD_LENGTH 8      /* word length is 8 bytes */
#define MAXOPFILES 100
static int ofp[MAXOPFILES];

int zop_(int *ioc, char *name, int *nsize, int *idisk,
         int *icode, int *ilab, int *name_length)
{
  char name_wrk[81];
  int length;
  int i, ibeg;
  char *c;

/* Find the first non-blank character in string "name" blank=32*/
  i = 0;
  c = (name+i);
  while ((int)*c == 32 || (int)*c == 0) {
     /* printf("i=%d   c=%c**  c=%d\n",i,*c,*c); */
     i++;
     c = (name+i);
  }
  ibeg = i;
  while ((int)*c != 32 && (int)*c != 0 && i < 80) {
     /* printf("i=%d   c=%c**  c=%d\n",i,*c,*c); */
     i++;
     c = (name+i);
  }
  length = i;

  strncpy(name_wrk, name, i);
  strncpy((name_wrk+i), "\0", 1);

/*  printf("%s   length=%d\n", name_wrk, name_length); */

  ofp[*ioc] = open(name_wrk, O_RDWR | O_CREAT, 0666);
  printf("file descriptor = %d\n\n", ofp[*ioc]);

  return 0;
}
   
int zcl_(int *ioc)
{
  close(ofp[*ioc]);
  return 0;
}

int zwr_(int *ioc, char *a, int *nwords, int *adres, int *lgivup, int *irr)
{
  int num_write, ncurr, offsets, nbytes;

  nbytes = *nwords * WORD_LENGTH;
/*  ncurr = ftell(ofp[*ioc]); */
/*  offsets = *adres * WORD_LENGTH; */
  offsets = (*adres - 0) * WORD_LENGTH;
  lseek(ofp[*ioc], offsets, SEEK_SET);
/*  num_write = fwrite(a, sizeof(char), *nwords, ofp[*ioc]); */
  num_write = write(ofp[*ioc], a, nbytes);
/*  printf("nwords=%d     number of words written=%d\n", *nwords, num_write);
*/

  return 0;
}

int zrd_(int *ioc, char *a, int *nwords, int *adres, int *lgivup, int *irr)
{
  int num_read, ncurr, offsets, nbytes;

  nbytes = *nwords * WORD_LENGTH;
/*  ncurr = ftell(ofp[*ioc]); */  
  offsets = (*adres - 0 ) * WORD_LENGTH;
  /*  printf("In zrdc: offsets = %d    ioc=%d   ncurr=%d\n", offsets, *ioc, ncurr);
*/
  lseek(ofp[*ioc], offsets, SEEK_SET);
/*  num_read = read(a, sizeof(char), *nwords, ofp[*ioc]); */
  num_read = read(ofp[*ioc], a, nbytes);
/*  printf("nwords=%d     number of words read=%d\n", *nwords, num_read);
*/

  return 0;
}

