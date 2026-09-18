/* print the C standard the compiler claims to support */

#include <stdio.h>
#include <stdlib.h>

int main(void)
{

printf("__STDC_VERSION__ %ld\n", __STDC_VERSION__);

return EXIT_SUCCESS;
}
