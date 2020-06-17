#!/usr/bin/env python3
"""
this script is because some systems may not see the Python entry_point via CMake.
This seemed like a quick workaround / backup for those cases.
"""

import gemini3d.compare

gemini3d.compare.cli()
