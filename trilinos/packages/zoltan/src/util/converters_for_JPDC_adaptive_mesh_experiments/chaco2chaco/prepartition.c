

/* ==== Put this code in your main timestep (=ts) loop AFTER you read in the new mesh ===*/

/*=== preFileOpen is an int ´flag´ ===*/

    /*=== preFileName ===*/
    sprintf(preFileName,"%s.pre",pio_info.pexo_fname);
    
    if (ts>0 && Proc==0)
      {
	preFileOpen=1;
	if ((preFile=fopen(preFileName,"r"))==NULL) preFileOpen=0;
	if (preFileOpen) 
	  {
	    readOK=readPrePartition(mesh, oldMesh, preFile);
	  }
	if ((preFileOpen==0) || (readOK==0))
	  {
	    prePartition(mesh, oldMesh, preFileName);
	  }
      }
/*=========================================================================================*/






/*=== Below are the functions called in the code above. Put it where you see fit:-) ===*/

#define SQUARE(x) ((x)*(x))


#define SQUAREDDISTANCE(DIMS, COORDS1, COORDS2, RES){	\
  d1=COORDS1[0]-COORDS2[0];\
  RES=SQUARE(d1);\
  if (DIMS>1)\
    {\
      d1=COORDS1[1]-COORDS2[1];\
      RES+=SQUARE(d1);\
      if (DIMS>2)\
	{\
	  d1=COORDS1[2]-COORDS2[2];\
	  RES+=SQUARE(d1);\
	}\
    }\
}



static void prePartition(MESH_INFO_PTR mesh, MESH_INFO_PTR oldMesh, char *fileName)
{

  int i,j, oldElems, newElems, closestIndex;
  double distSq, closestSoFar;
  double d1,d2;
  FILE *preFile;

  /*=== For each element in mesh, find the closest element in oldMesh ===*/
  /*=== Assign part...===*/

  if ((preFile=fopen(fileName,"w"))==NULL) {printf("Cannot open pre partitioning file\n");return;}
  
  newElems=mesh->num_elems;oldElems=oldMesh->num_elems;
  for (i=0; i<newElems; i++)
    {
      closestSoFar=100000.0;closestIndex=0;
      for (j=0; j<oldElems; j++)
	{
	  SQUAREDDISTANCE(mesh->num_dims, mesh->elements[i].avg_coord, oldMesh->elements[j].avg_coord, distSq);
	  if (distSq<closestSoFar)
	    {
	      closestSoFar=distSq;closestIndex=j;
	    }
	}
      mesh->elements[i].my_part=oldMesh->elements[closestIndex].my_part;
      fprintf(preFile,"%d\n",closestIndex);
    }
  fclose(preFile);
}


static int readPrePartition(MESH_INFO_PTR mesh, MESH_INFO_PTR oldMesh, FILE *preFile)
{

  int j, newElems, closestIndex, num;
  char line[40];

  /*=== For each element in mesh, read the closest element in oldMesh ===*/

  newElems=mesh->num_elems;
  for (j=0; j<newElems; j++)
    {
       if (!(fgets(line, 40, preFile))) 
	{
	  printf("Could not read the line with DID and nadj from th prepartitioning file\n");
	  return 0;
	}
      num=sscanf(line,"%d",&closestIndex);
      if (num!=1) {printf("Could not read closest from the pre partitioning file!\n"); return 0;}
      if ( (closestIndex > mesh->num_elems) && (closestIndex<0) )
	{printf("Baaad closest index (%d)\n",closestIndex);return 0;}
      mesh->elements[j].my_part=oldMesh->elements[closestIndex].my_part;
    }
  fclose(preFile);
  return 1;
}
