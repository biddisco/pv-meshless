/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkTipsyReader.cxx,v $
=========================================================================*/
#include "vtkTipsyReader.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h" 
#include "vtkIntArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDistributedDataFilter.h"
#include "vtkMultiProcessController.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkDataArraySelection.h"
#include <cmath>
#include <assert.h>

#include "FileSeriesFinder.h"

vtkCxxRevisionMacro(vtkTipsyReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkTipsyReader);
//----------------------------------------------------------------------------
vtkSmartPointer<vtkFloatArray> AllocateDataArray(
  vtkDataSet *output, const char* arrayName, int numComponents, unsigned long numTuples)
{
  vtkSmartPointer<vtkFloatArray> dataArray=vtkSmartPointer<vtkFloatArray>::New();
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetNumberOfTuples(numTuples);
	dataArray->SetName(arrayName);
	// initializes everything to zero
	for (int i=0; i < numComponents; ++i) {
		dataArray->FillComponent(i, 0.0);
  }
  output->GetPointData()->AddArray(dataArray);
  return dataArray;
}

//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
  this->MarkFileName      = 0; // this file is optional
  this->FileName          = 0;
	this->DistributeDataOn  = 1;
  this->UpdatePiece       = 0;
  this->UpdateNumPieces   = 0;
  this->SetNumberOfInputPorts(0); 
  //
  this->PointDataArraySelection  = vtkDataArraySelection::New();
  //
  this->Positions     = NULL;
  this->Vertices      = NULL;
  this->GlobalIds     = NULL;
  this->ParticleIndex = 0;
  this->Potential   = NULL;
  this->Mass        = NULL;
  this->EPS         = NULL;
  this->RHO         = NULL;
  this->Hsmooth     = NULL;
  this->Temperature = NULL;
  this->Metals      = NULL;
  this->Tform       = NULL;
  this->Velocity    = NULL;
  //
  this->Finder      = NULL;
}

//----------------------------------------------------------------------------
vtkTipsyReader::~vtkTipsyReader()
{
  this->SetFileName(0);
  this->SetMarkFileName(0);
	this->SetDistributeDataOn(0);
  this->PointDataArraySelection->Delete();
}

//----------------------------------------------------------------------------
void vtkTipsyReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n"
		 << indent << "MarkFileName: "
		 << (this->MarkFileName ? this->MarkFileName : "(none)") << "\n";
}

//----------------------------------------------------------------------------
void vtkTipsyReader::GetTimeStepValues(std::vector<double> &steps)
{
  steps = this->TimeStepValues;
}

//----------------------------------------------------------------------------
TipsyHeader vtkTipsyReader::ReadTipsyHeader(ifTipsy& tipsyInfile)
{
	// Initializing for reading  
	TipsyHeader       h; 
	// Reading in the header
  tipsyInfile >> h;
	return h;
}

//----------------------------------------------------------------------------
vtkstd::vector<unsigned long> vtkTipsyReader::ReadMarkedParticleIndices(
	TipsyHeader& tipsyHeader,ifTipsy& tipsyInfile)
{
	ifstream markInFile(this->MarkFileName);
	vtkstd::vector<unsigned long> markedParticleIndices;
	if(!markInFile)
 		{
 		vtkErrorMacro("Error opening marked particle file: " 
									<< this->MarkFileName 
									<< " please specify a valid mark file or none at all.\
									 For now reading all particles.");
		return markedParticleIndices;
 		}
	else
		{
		unsigned long mfIndex,mfBodies,mfGas,mfStar,mfDark,numBodies;
		// first line of the mark file is of a different format:
		// intNumBodies intNumGas intNumStars
		if(markInFile >> mfBodies >> mfGas >> mfStar)
			{
	 		mfDark=mfBodies-mfGas-mfStar;
			if(mfBodies!=tipsyHeader.h_nBodies || mfDark!=tipsyHeader.h_nDark \
				|| mfGas!=tipsyHeader.h_nSph || mfStar!=tipsyHeader.h_nStar)
	 			{
	 			vtkErrorMacro("Error opening marked particle file, wrong format,\
	 										number of particles do not match Tipsy file: " 
											<< this->MarkFileName 
											<< " please specify a valid mark file or none at all.\
											For now reading all particles.");
	 			}
			else
				{
				// The rest of the file is is format: intMarkIndex\n
				// Read in the rest of the file file, noting the marked particles
				while(markInFile >> mfIndex)
					{
					// Insert the next marked particle index into the vector
					// subtracting one as the indices in marked particle file
					// begin at 1, not 0
					markedParticleIndices.push_back(mfIndex-1);
					}
				// closing file
				markInFile.close();
				// read file successfully
				// vtkDebugMacro("Read " << numBodies<< " marked point indices.");
				}	
	 		}
		}
		//empty if none were read, otherwise size gives number of marked particles
		return markedParticleIndices;
}

//----------------------------------------------------------------------------
void vtkTipsyReader::ReadAllParticles(TipsyHeader& tipsyHeader,
	ifTipsy& tipsyInfile,int piece,int numpieces,vtkPolyData* output)
{
	unsigned long pieceSize = floor(tipsyHeader.h_nBodies*1./numpieces);
	unsigned long beginIndex = piece*pieceSize;
	unsigned long endIndex = (piece == numpieces - 1) ? \
	 	tipsyHeader.h_nBodies : (piece+1)*pieceSize;
	// Allocates vtk scalars and vector arrays to hold particle data, 
	this->AllocateAllTipsyVariableArrays(endIndex-beginIndex,output);
	for(unsigned long i=beginIndex; i < endIndex; i++)  
  	{ 
		this->ReadParticle(i,tipsyHeader,tipsyInfile,output);
  	}
}

//----------------------------------------------------------------------------
void vtkTipsyReader::ReadMarkedParticles(
	vtkstd::vector<unsigned long>& markedParticleIndices,
	TipsyHeader& tipsyHeader,
	ifTipsy& tipsyInfile,
	vtkPolyData* output)
{
	// Allocates vtk scalars and vector arrays to hold particle data, 
	// As marked file was read, only allocates numBodies which 
	// now equals the number of marked particles
	this->AllocateAllTipsyVariableArrays(markedParticleIndices.size(),output);
	unsigned long nextMarkedParticleIndex=0;
	for(vtkstd::vector<unsigned long>::iterator it = \
	 	markedParticleIndices.begin(); it != markedParticleIndices.end(); ++it)		
		{
 		nextMarkedParticleIndex=*it;
		// reading in the particle
		this->ReadParticle(nextMarkedParticleIndex,tipsyHeader,
			tipsyInfile,output);
		}
}

//----------------------------------------------------------------------------
tipsypos::section_type vtkTipsyReader::SeekToIndex(unsigned long index,
	TipsyHeader& tipsyHeader, ifTipsy& tipsyInfile)
{
	if(index < tipsyHeader.h_nSph)
		{
		// we are seeking a gas particle
		tipsyInfile.seekg(tipsypos(tipsypos::gas,index));	
		return tipsypos::gas;
		}
	else if(index < tipsyHeader.h_nDark+tipsyHeader.h_nSph)
		{
		// we are seeking a dark particle, which begin at index 0 after
		// gas particles
		index-=tipsyHeader.h_nSph;
		tipsyInfile.seekg(tipsypos(tipsypos::dark,index));
		return tipsypos::dark;
		}
	else if(index < tipsyHeader.h_nBodies)
		{
		// we are seeking a star particle, which begin at index zero after gas
		// and dark particles
		index-=(tipsyHeader.h_nSph+tipsyHeader.h_nDark);
		tipsyInfile.seekg(tipsypos(tipsypos::star,index));	
		return tipsypos::star;
		}
	else
		{
		vtkErrorMacro("An index is greater than the number of particles in the 	file, unable to read");
		return tipsypos::invalid;
		}
}

//----------------------------------------------------------------------------
//i,tipsyHeader,tipsyInfile,output)
vtkIdType vtkTipsyReader::ReadParticle(unsigned long index, 
	TipsyHeader& tipsyHeader, ifTipsy& tipsyInfile, vtkPolyData* output) 
{
  // allocating variables for reading
  vtkIdType id;
  TipsyGasParticle  g;
  TipsyDarkParticle d;
  TipsyStarParticle s;
	tipsypos::section_type particleSection=this->SeekToIndex(index,
		tipsyHeader,tipsyInfile);
  static bool error_reported = false;
	switch(particleSection)
		{
   case tipsypos::gas:
     tipsyInfile >> g;
		 id=this->ReadGasParticle(output,g);
     break;
   case tipsypos::dark:
     tipsyInfile >> d;
		 id=this->ReadDarkParticle(output,d);
     break;
   case tipsypos::star:
     tipsyInfile >> s;
		 id=this->ReadStarParticle(output,s);
     break;
   default:
     // only do this once
     if (!error_reported) {
      vtkErrorMacro(<<"Wrong particle type - not supported")
      error_reported = true;
     }
     break;
    }
  this->GlobalIds->SetValue(id,index);
  return id;
}

//----------------------------------------------------------------------------
vtkIdType vtkTipsyReader::ReadBaseParticle(vtkPolyData* output,
	TipsyBaseParticle& b) 
{
  this->Positions->SetPoint(ParticleIndex, b.pos);
  //
  if (this->Velocity)  this->Velocity->SetTuple(ParticleIndex, b.vel);
  if (this->Mass)      this->Mass->SetTuple1(ParticleIndex, b.mass);
  if (this->Potential) this->Potential->SetTuple1(ParticleIndex, b.phi);
	return ParticleIndex++;
}

//----------------------------------------------------------------------------
vtkIdType vtkTipsyReader::ReadGasParticle(vtkPolyData* output,
 	TipsyGasParticle& g)
{
 	vtkIdType id=ReadBaseParticle(output,g);
  if (this->RHO)         this->RHO->SetTuple(id, &g.rho);
  if (this->Temperature) this->Temperature->SetTuple1(id, g.temp);
  if (this->Hsmooth)     this->Hsmooth->SetTuple1(id, g.hsmooth);
  if (this->Metals)      this->Metals->SetTuple1(id, g.metals);
	if (this->Type)      this->Type->SetTuple1(id, 0.0);
	return id;
}

//----------------------------------------------------------------------------
vtkIdType vtkTipsyReader::ReadStarParticle(vtkPolyData* output,
 	TipsyStarParticle& s)
{
 	vtkIdType id=ReadBaseParticle(output,s);
  if (this->EPS)    this->EPS->SetTuple1(id, s.eps);
  if (this->Metals) this->Metals->SetTuple1(id, s.metals);
  if (this->Tform)  this->Tform->SetTuple1(id, s.tform);
	if (this->Type)      this->Type->SetTuple1(id, 2.0);
	return id;
}

//----------------------------------------------------------------------------	
vtkIdType vtkTipsyReader::ReadDarkParticle(vtkPolyData* output,
	TipsyDarkParticle& d)
{
 	vtkIdType id=this->ReadBaseParticle(output,d);
  if (this->EPS) this->EPS->SetTuple1(id, d.eps);
	if (this->Type)      this->Type->SetTuple1(id, 1.0);
	return id;
}
		
//----------------------------------------------------------------------------
void vtkTipsyReader::AllocateAllTipsyVariableArrays(vtkIdType numBodies,
	vtkPolyData* output)
{
  // Allocate objects to hold points and vertex cells. 
  this->Positions = vtkSmartPointer<vtkPoints>::New();
  this->Positions->SetDataTypeToFloat();
  this->Positions->SetNumberOfPoints(numBodies);
  //
  this->Vertices  = vtkSmartPointer<vtkCellArray>::New();
  vtkIdType *cells = this->Vertices->WritePointer(numBodies, numBodies*2);
  for (vtkIdType i=0; i<numBodies; ++i) {
    cells[i*2]   = 1;
    cells[i*2+1] = i;
  }

  //
  this->GlobalIds = vtkSmartPointer<vtkIdTypeArray>::New();
  this->GlobalIds->SetName("global_id");
  this->GlobalIds->SetNumberOfTuples(numBodies);

	// Storing the points and cells in the output data object.
  output->SetPoints(this->Positions);
  output->SetVerts(this->Vertices); 

  // allocate velocity first as it uses the most memory and on my win32 machine 
  // this helps load really big data without alloc failures.
  if (this->GetPointArrayStatus("Velocity")) 
    this->Velocity = AllocateDataArray(output,"velocity",3,numBodies);
  else 
    this->Velocity = NULL;
  if (this->GetPointArrayStatus("Potential")) 
    this->Potential = AllocateDataArray(output,"potential",1,numBodies);
  else 
    this->Potential = NULL;
  if (this->GetPointArrayStatus("Mass"))
    this->Mass = AllocateDataArray(output,"mass",1,numBodies);
  else 
    this->Mass = NULL;
  if (this->GetPointArrayStatus("Eps")) 
    this->EPS = AllocateDataArray(output,"eps",1,numBodies);
  else 
    this->EPS = NULL;
  if (this->GetPointArrayStatus("Rho")) 
    this->RHO = AllocateDataArray(output,"rho",1,numBodies);
  else 
    this->RHO = NULL;
  if (this->GetPointArrayStatus("Hsmooth")) 
    this->Hsmooth = AllocateDataArray(output,"hsmooth",1,numBodies);
  else 
    this->Hsmooth = NULL;
  if (this->GetPointArrayStatus("Temperature"))
    this->Temperature = AllocateDataArray(output,"temperature",1,numBodies);
  else 
    this->Temperature = NULL;

  if (this->GetPointArrayStatus("Metals"))
    this->Metals = AllocateDataArray(output,"metals",1,numBodies);
  else 
    this->Metals = NULL;
  if (this->GetPointArrayStatus("Tform"))
    this->Tform = AllocateDataArray(output,"tform",1,numBodies);
  else 
    this->Tform = NULL;
	if (this->GetPointArrayStatus("Type"))
    this->Type = AllocateDataArray(output,"type",1,numBodies);
  else 
    this->Type = NULL;
}
//----------------------------------------------------------------------------
int vtkTipsyReader::RequestInformation(
	vtkInformation* vtkNotUsed(request),
	vtkInformationVector** vtkNotUsed(inputVector),
	vtkInformationVector* outputVector)
{
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	// means that the data set can be divided into an arbitrary number of pieces
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
		-1);

  this->PointDataArraySelection->AddArray("Potential");
  this->PointDataArraySelection->AddArray("Mass");
  this->PointDataArraySelection->AddArray("Eps");
  this->PointDataArraySelection->AddArray("Rho");
  this->PointDataArraySelection->AddArray("Hsmooth");
  this->PointDataArraySelection->AddArray("Temperature");
  this->PointDataArraySelection->AddArray("Metals");
  this->PointDataArraySelection->AddArray("Tform");
	this->PointDataArraySelection->AddArray("Type");
  this->PointDataArraySelection->AddArray("Velocity");

  if (!this->Finder) {   
    this->Finder = new FileSeriesFinder();
    this->Finder->SetPrefixRegEx("(.*\\.)");
    this->Finder->SetTimeRegEx("([0-9]+)");
    this->Finder->Scan(this->FileName);
    this->Finder->GetTimeValues(this->TimeStepValues);
    this->NumberOfTimeSteps = this->Finder->GetNumberOfTimeSteps();
  }


	return 1;
}
/*
* Reads a file, optionally only the marked particles from the file, 
* in the following order:
* 1. Open Tipsy binary
* 2. Read Tipsy header (tells us the number of particles of each type we are 
*    dealing with)
* NOTE: steps 3, 5 are currently not parallel
* 3. Read mark file indices from marked particle file, if there is one
* 4. Read either marked particles only or all particles
* 5. If an attribute file is additionally specified, reads this additional
* 	 attribute into a data array, reading only those marked if necessary.
*/
//----------------------------------------------------------------------------
int vtkTipsyReader::RequestData(vtkInformation*,
	vtkInformationVector**,vtkInformationVector* outputVector)
{
  //
	// Make sure we have a file to read.
  //
  if(!this->FileName)
	  {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
    }

	// Open the tipsy standard file and abort if there is an error.
	ifTipsy tipsyInfile;
  tipsyInfile.open(this->FileName,"standard");
  if (!tipsyInfile.is_open()) 
		{
	  vtkErrorMacro("Error opening file " << this->FileName);
	  return 0;	
    }

  // Get output information
	vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the output polydata
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkSmartPointer<vtkPolyData> tipsyReadInitialOutput = vtkSmartPointer<vtkPolyData>::New();

  // get this->UpdatePiece information
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
	this->UpdateNumPieces =outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  // reset counter before reading
  this->ParticleIndex = 0;

  // Read the header from the input
	TipsyHeader tipsyHeader=this->ReadTipsyHeader(tipsyInfile);
	// Next considering whether to read in a mark file, 
	// and if so whether that reading was a success 
	vtkstd::vector<unsigned long> markedParticleIndices;
	if(this->MarkFileName && strcmp(this->MarkFileName,"")!=0)
		{
		// Reading only marked particles
		// Make sure we are not running in parallel, this filter does not work in 
		// parallel
		if(this->UpdateNumPieces>1)
			{
			vtkErrorMacro("Reading from a mark file is not supported in parallel.");
			return 0;
			}
		vtkDebugMacro("Reading marked point indices from file:" 
			<< this->MarkFileName);
		markedParticleIndices=this->ReadMarkedParticleIndices(tipsyHeader,
			tipsyInfile);
		}
  // Read every particle and add their position to be displayed, 
	// as well as relevant scalars
	if(markedParticleIndices.empty())
		{
		// no marked particle file or there was an error reading the mark file, 
		// so reading all particles
		vtkDebugMacro("Reading all points from file " << this->FileName);
		this->ReadAllParticles(tipsyHeader,tipsyInfile, this->UpdatePiece, this->UpdateNumPieces, tipsyReadInitialOutput);
		}
	else 
		{
		//reading only marked particles
		assert(this->UpdateNumPieces==1);
		vtkDebugMacro("Reading only the marked points in file: " \
				<< this->MarkFileName << " from file " << this->FileName);
		this->ReadMarkedParticles(markedParticleIndices, tipsyHeader,tipsyInfile, tipsyReadInitialOutput);	
		}
  // Close the tipsy in file.
	tipsyInfile.close();
	// If we need to, run D3 on the tipsyReadInitialOutput
	// producing one level of ghost cells
	if (this->GetDistributeDataOn() && this->UpdateNumPieces>1)
		{
		vtkSmartPointer<vtkDistributedDataFilter> d3 = \
		    vtkSmartPointer<vtkDistributedDataFilter>::New();
		d3->AddInput(tipsyReadInitialOutput);
		d3->UpdateInformation();
		d3->Update();
		// Changing output to output of d3
	 	output->ShallowCopy(d3->GetOutput()); 

    // create new vertices because the D3 outputs UnstructuredGrid
    vtkIdType N = output->GetNumberOfPoints();
    this->Vertices  = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType *cells = this->Vertices->WritePointer(N, N*2);
    for (vtkIdType i=0; i<N; ++i) {
      cells[i*2]   = 1;
      cells[i*2+1] = i;
    }
    output->SetVerts(this->Vertices);
		}
	else
		{
		output->ShallowCopy(tipsyReadInitialOutput);
		}
	// Read Successfully
	vtkDebugMacro("Read " << output->GetPoints()->GetNumberOfPoints() \
		<< " points.");
  //
  // release memory smartpointers - just to play safe.
  this->Vertices    = NULL;
  this->GlobalIds   = NULL;
  this->Positions   = NULL;
  this->Potential   = NULL;
  this->Mass        = NULL;
  this->EPS         = NULL;
  this->RHO         = NULL;
  this->Hsmooth     = NULL;
  this->Temperature = NULL;
  this->Metals      = NULL;
  this->Tform       = NULL;
	this->Type       = NULL;
  this->Velocity    = NULL;
  //
 	return 1;
}
//----------------------------------------------------------------------------
// Below : Boiler plate code to handle selection of point arrays
//----------------------------------------------------------------------------
const char* vtkTipsyReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkTipsyReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkTipsyReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name)) {
    if (status) {
      this->PointDataArraySelection->EnableArray(name);
    }
    else {
      this->PointDataArraySelection->DisableArray(name);
    }
    this->Modified();
  }
}
//----------------------------------------------------------------------------
void vtkTipsyReader::Enable(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}
//----------------------------------------------------------------------------
void vtkTipsyReader::Disable(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}
//----------------------------------------------------------------------------
void vtkTipsyReader::EnableAll()
{
  this->PointDataArraySelection->EnableAllArrays();
}
//----------------------------------------------------------------------------
void vtkTipsyReader::DisableAll()
{
  this->PointDataArraySelection->DisableAllArrays();
}
//----------------------------------------------------------------------------
int vtkTipsyReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
