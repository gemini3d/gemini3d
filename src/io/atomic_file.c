/* SPDX-License-Identifier: Apache-2.0 */
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <string.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

int gemini_stage_create(char *path) {
#ifdef _WIN32
    if (_mktemp_s(path, strlen(path) + 1)) return -1;
    return _open(path, _O_CREAT | _O_EXCL | _O_RDWR | _O_BINARY, _S_IREAD | _S_IWRITE);
#else
    return mkstemp(path);
#endif
}
int gemini_stage_close(int fd) {
#ifdef _WIN32
    return _close(fd);
#else
    return close(fd);
#endif
}
int gemini_stage_publish(const char *old_path, const char *new_path) {
#ifdef _WIN32
    return MoveFileExA(old_path, new_path, MOVEFILE_REPLACE_EXISTING) ? 0 : -1;
#else
    return rename(old_path, new_path);
#endif
}
