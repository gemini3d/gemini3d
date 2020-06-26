---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

If it's a problem getting Gemini to build (compile), please let us know the output of these commands (stopping at the command that fails)

```sh
cmake -B build
```

```sh
cmake --build build
```

```sh
cd build

ctest -V
```
