#if __STDC_VERSION__ < 202311L
#include <stdbool.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ffilesystem.h"

#ifdef _MSC_VER
#include <crtdbg.h>
#endif


int main(void){

#ifdef _MSC_VER
  _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
  _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
  _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);
#endif

  size_t N = 5;
  char* buf = (char*) malloc(N);
  int i=0;
  size_t L;

  if(fs_normal("abcedf", buf, N) != 0){
    fprintf(stderr, "ERROR: fs_normal(abcdef) did not handle overflow properly\n");
    fprintf(stderr, "       buf = %s\n", buf);
    i++;
  }

  if(fs_expanduser("~", buf, N) != 0){
    fprintf(stderr, "ERROR: fs_expanduser(~) did not handle overflow properly\n");
    fprintf(stderr, "       buf = %s\n", buf);
    i++;
  }

  if(fs_expanduser("abcedf", buf, N) != 0){
    fprintf(stderr, "ERROR: fs_expanduser(abcdef) did not handle overflow properly\n");
    fprintf(stderr, "       buf = %s\n", buf);
    i++;
  }

  L = fs_file_name("abcedf", buf, N);
  if(L != 0){
    fprintf(stderr, "ERROR: fs_file_name(abcdef) did not handle overflow properly\n");
    fprintf(stderr, "L = %zu       buf = %s\n", L, buf);
    i++;
  }

  L = fs_suffix("abcedf", buf, N);
  if(L != 0){
    fprintf(stderr, "ERROR: fs_stem(abcdef) did not handle overflow properly\n");
    fprintf(stderr, "L = %zu       buf = %s\n", L, buf);
    i++;
  }
  L = fs_parent("abcedf/a", buf, N);
  if(L != 0){
    fprintf(stderr, "ERROR: fs_parent(abcdef) did not handle overflow properly\n");
    fprintf(stderr, "L = %zu       buf = %s\n", L, buf);
    i++;
  }

  if(fs_with_suffix(".abcedf", "txt", buf, N) != 0){
    fprintf(stderr, "ERROR: fs_with_suffix(abcdef, txt) did not handle overflow properly\n");
    fprintf(stderr, "       buf = %s\n", buf);
    i++;
  }

  if(fs_root("abcdef", buf, N) != 0){
    fprintf(stderr, "ERROR: fs_root(abcdef) did not handle overflow properly\n");
    fprintf(stderr, "       buf = %s\n", buf);
    i++;
  }

  if(fs_absolute("abcdef", "zyxwvu", buf, N) != 0){
    fprintf(stderr, "ERROR: fs_absolute(abcdef) did not handle overflow properly\n");
    fprintf(stderr, "       buf = %s\n", buf);
    i++;
  }

  if(fs_longname("abcdef", buf, N) != 0){
    fprintf(stderr, "ERROR: fs_longname(abcdef) did not handle overflow properly\n");
    fprintf(stderr, "       buf = %s\n", buf);
    i++;
  }

  if(fs_shortname("abcdef", buf, N) != 0){
    fprintf(stderr, "ERROR: fs_shortname(abcdef) did not handle overflow properly\n");
    fprintf(stderr, "       buf = %s\n", buf);
    i++;
  }

  free(buf);

  if(i > 0)
    return EXIT_FAILURE;

  printf("OK: check overflow behavior\n");

  return EXIT_SUCCESS;
}
