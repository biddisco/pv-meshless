/*=========================================================================

  Project: meshless plugin

  Copyright (c) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing this
  copyright notice appears on all copies of source code and an
  acknowledgment appears with any substantial usage of the code.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  This code is derived from an earlier work and is distributed
  with permission from, and thanks to

=========================================================================*/
/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkTipsyReader.cxx,v $
=========================================================================*/

/*=========================================================================

  Program:   VTK/ParaView Los Alamos National Laboratory Modules (PVLANL)
  Module:    $RCSfile: vtkPTipsyReader.h,v $

Copyright (c) 2009 Los Alamos National Security, LLC

All rights reserved.

Copyright 2009. Los Alamos National Security, LLC. 
This software was produced under U.S. Government contract DE-AC52-06NA25396 
for Los Alamos National Laboratory (LANL), which is operated by 
Los Alamos National Security, LLC for the U.S. Department of Energy. 
The U.S. Government has rights to use, reproduce, and distribute this software. 
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,
EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
If software is modified to produce derivative works, such modified software 
should be clearly marked, so as not to confuse it with the version available 
from LANL.
 
Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
-   Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer. 
-   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution. 
-   Neither the name of Los Alamos National Security, LLC, Los Alamos National
    Laboratory, LANL, the U.S. Government, nor the names of its contributors
    may be used to endorse or promote products derived from this software 
    without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
// .NAME vtkPTipsyReader - Read a binary Tipsy cosmology data file
//
// .SECTION Description
// vtkPTipsyReader creates a vtkUnstructuredGrid from a binary cosmology file.
// The file contains fields for:
//     x_position, x_velocity (float)
//     y_position, y_velocity (float)
//     z-position, z_velocity (float)
//     mass (float)
//     identification tag (integer)
//
// If the file contains particle information x,y,z is the location of the
// particle in simulation space with a velocity vector and a mass which
// will be the same for all particles.
//

#ifndef __vtkPTipsyReader_h
#define __vtkPTipsyReader_h

#include "vtkPCosmoReader.h"

class VTK_EXPORT vtkPTipsyReader : public vtkPCosmoReader
{
public:
  static vtkPTipsyReader *New();
  vtkTypeMacro(vtkPTipsyReader, vtkPCosmoReader);

  // Description:
  // Get/Set whether to distribute data
	vtkSetMacro(RecentreBoundingBox,int);
	vtkGetMacro(RecentreBoundingBox,int);
	vtkBooleanMacro(RecentreBoundingBox,int);

protected:
   vtkPTipsyReader();
  ~vtkPTipsyReader();
 
  virtual int RequestInformation
    (vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestData
    (vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  int RecentreBoundingBox;

private:
  vtkPTipsyReader(const vtkPTipsyReader&);  // Not implemented.
  void operator=(const vtkPTipsyReader&);  // Not implemented.
};

#endif

