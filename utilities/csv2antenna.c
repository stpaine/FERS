// Create a FERS antenna description file from 2 sets of CSV antenna data
// Marc Brooker mbrooker@rrsg.ee.uct.ac.za
// 12 June 2007

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ProcessCSV(FILE *fpin, FILE *fpout, const char *enc_tag, const char *d1_tag, const char *d2_tag)
{
  char *buff = malloc(2048);
  char *line;
  line = fgets(buff, 2048, fpin);
  while (line != NULL) {
    char *comma = strchr(line, ',');
    if (!comma) {
      fprintf(stderr, "[ERROR] Malformed CSV line in input file %s\n\n", line);
      exit(1);
    }
    *comma = 0;
    comma++;
    comma[strlen(comma)-1]=0;
    fprintf(stderr, "%s\n", line);
    fprintf(stderr, "%s\n", comma+1);
    fprintf(fpout, "\t<%s>\n\t\t<%s>%s</%s><%s>%s</%s>\n\t</%s>\n", enc_tag, d1_tag, line, d1_tag, d2_tag, comma, d2_tag, enc_tag);
    line = fgets(line, 2048, fpin);
  }
}

int main(int argc, char *argv[])
{
  FILE *fpout, *fpin1, *fpin2;
  if (argc != 4) {
    fprintf(stderr, "Usage: csv2antenna <outfile> <elevation gains> <azimuth gains>\n");
    exit(2);
  }
  fpin1 = fopen(argv[2], "r");
  if (!fpin1) {
    perror("Could not open input file");
    exit(2);
  }
  fpin2 = fopen(argv[3], "r");
  if (!fpin2) {
    perror("Could not open input file");
    exit(2);
  }
  fpout = fopen(argv[1], "w");
  if (!fpout) {
    perror("Could not open input file");
    exit(2);
  }
  fprintf(fpout, "<antenna>\n<elevation>\n");
  ProcessCSV(fpin1, fpout, "gainsample", "angle", "gain");
  fprintf(fpout, "</elevation>\n<azimuth>\n");
  ProcessCSV(fpin2, fpout, "gainsample", "angle", "gain");
  fprintf(fpout, "</azimuth>\n</antenna>\n");
  fclose(fpin1);
  fclose(fpin2);
  fclose(fpout);
}

    
  
    
  
