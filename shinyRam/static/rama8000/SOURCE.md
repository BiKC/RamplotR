# Rama8000 reference tables

These six 180 x 180 score tables are vendored from the cctbx/Phenix Rama8000 implementation used by `mmtbx.validation.ramalyze`.

Source repository: `cctbx/cctbx_project`
Source commit: `b897de9a4d9f0a81945a05158c915ff8b3d37b4a`
Generated header: `mmtbx/validation/ramachandran/rama8000_tables.h`
Header blob SHA: `b1e8812b58b802b73333dd3b391771c22e08bc82`

Residue classes:
- general
- glycine
- cis-proline
- trans-proline
- pre-proline
- isoleucine or valine

The tables are sampled at odd-numbered two-degree grid coordinates from -179 to 179 degrees. RamplotR uses periodic bilinear interpolation and the same score thresholds as cctbx `rama_eval.h`.

cctbx is distributed under a permissive BSD-style license. The relevant license and copyright notices are reproduced below.

## License notice

cctbx Copyright (c) 2006 - 2026, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy). All
rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the University of California, Lawrence Berkeley
   National Laboratory, U.S. Dept. of Energy nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THE SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
WITHOUT EXPRESS OR IMPLIED WARRANTIES. See the upstream `LICENSE.txt` for
the complete disclaimer.
