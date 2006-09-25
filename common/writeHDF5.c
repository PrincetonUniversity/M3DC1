#include <stdio.h>
#include <hdf5.h>

static char filename[] = "C1.h5";

/* Utilities modified from the AVS library */
static void write_int_attr( hid_t parent_id, const char *name, int value )
{
    hid_t h5_dspace_id, h5_attr_id;

    h5_dspace_id = H5Screate( H5S_SCALAR );
    if( h5_dspace_id >= 0 ) {
      /* First, see if this attribute already exists */
      h5_attr_id = H5Aopen_name(parent_id, name);
      if (h5_attr_id < 0) /* If not, create it */
        h5_attr_id = H5Acreate( parent_id, name, H5T_NATIVE_INT_g,
                                h5_dspace_id, H5P_DEFAULT );

        if( h5_attr_id >= 0 ) {
            H5Awrite( h5_attr_id, H5T_NATIVE_INT_g, &value );
            H5Aclose( h5_attr_id );	/* close attribute */
        }
        H5Sclose( h5_dspace_id );	/* close dataspace */
    }
    return;
}

static void write_double_attr(hid_t parent_id, const char *name, double value)
{
    hid_t h5_dspace_id, h5_attr_id;

    h5_dspace_id = H5Screate( H5S_SCALAR );
    if( h5_dspace_id >= 0 ) {
      /* First, see if this attribute already exists */
      h5_attr_id = H5Aopen_name(parent_id, name);
      if (h5_attr_id < 0) /* If not, create it */
        h5_attr_id = H5Acreate( parent_id, name, H5T_NATIVE_DOUBLE_g,
                                h5_dspace_id, H5P_DEFAULT );

        if( h5_attr_id >= 0 ) {
            H5Awrite(h5_attr_id, H5T_NATIVE_DOUBLE_g, &value);
            H5Aclose(h5_attr_id);	/* close attribute */
        }
        H5Sclose( h5_dspace_id );	/* close dataspace */
    }
    return;
}

static void write_string_attr( hid_t parent_id, const char *name,
                               const char *value )
{
    hid_t h5_dspace_id, h5_dtype_id, h5_attr_id;

    h5_dspace_id = H5Screate( H5S_SCALAR );
    if (h5_dspace_id > 0) {
        h5_dtype_id   = H5Tcopy( H5T_C_S1_g );
        if( h5_dtype_id > 0 ) {
            H5Tset_size( h5_dtype_id, strlen(value)+1 );
            h5_attr_id = H5Acreate ( parent_id, name,
                                     h5_dtype_id, h5_dspace_id, H5P_DEFAULT );
            if( h5_attr_id > 0 ) {
                H5Awrite( h5_attr_id, h5_dtype_id, (char *)value );
                H5Aclose( h5_attr_id );	/* close attribute */
            }
            H5Tclose( h5_dtype_id );	/* close datatype */
        }
        H5Sclose( h5_dspace_id );	/* close dataspace */
    }
    return;
}

/* My utilities (borrowing from AVS) */
int createhdf5_(double *xcoord, int *nx, double *ycoord, int *ny, int *nvar)
{
  int    fld_npoints = *nx + *ny; /* Points determining the mesh */
  int    fld_nspace = 2; /* Dimensionality of mesh */
  int    ii, fld_dims[2];
  double  *fld_points;
  herr_t status;
  hid_t  att_id, file_id, dataset, dataspace, root_id;
  hsize_t dims[2];
  H5E_auto_t olderrfunc;
  int i;


  /*printf("nx = %d, ny = %d, nvar = %d/n",*nx,*ny,*nvar);
  
     for (i=0;i<*nx;i++)  
    printf("%d %lf %lf\n",i,xcoord[i],ycoord[i]);*/


  /* Create a new HDF5 file, failing if such a file already exists */
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0) {
    fprintf(stderr, "Failed to create HDF5 file %s.\n", filename);
    return (int)file_id;
  }

  /* Temporarily turn off HDF5 error messages */
  H5Eget_auto(&olderrfunc, NULL);
  H5Eset_auto(NULL, NULL);

  /* Grab the root group, set some attributes for AVS to read */
  if ((root_id = H5Gopen(file_id, "/")) < 0) {
    H5Fclose(file_id);
    return (int)root_id;
  }
  /*write_string_attr(root_id, "XP_CLASS", "Mesh_Rect+Node_Data");*/
  write_string_attr(root_id, "XP_CLASS", "Mesh_Rect+Time_Node_Data");
  write_int_attr(root_id, "grid_type", 2 );
  write_int_attr(root_id, "npoints", fld_npoints);
  write_int_attr(root_id, "nspace", fld_nspace);
  write_int_attr(root_id, "ndim", fld_nspace);
  write_int_attr(root_id, "nnode_data", *nvar); /* vars per slice */

  /* Set dimensions */
  dims[0] = fld_nspace;
  dataspace = H5Screate_simple(1, dims, NULL);
  if (dataspace >= 0) {
    /* Create an attribute */
    att_id = H5Acreate(root_id, "dims", H5T_NATIVE_INT_g, dataspace,
		       H5P_DEFAULT);
    if (att_id >= 0) {
      fld_dims[0] = *nx;
      fld_dims[1] = *ny;
      H5Awrite(att_id, H5T_NATIVE_INT_g, fld_dims);
      H5Aclose(att_id);	/* close attribute */
    }
    H5Sclose(dataspace);
  }

  /* Set point coordinates */
  dims[0] = fld_npoints;
  dims[1] = fld_nspace;
  dataspace = H5Screate_simple(2, dims, NULL);
  if (dataspace >= 0) {
    /* Create the 2D dataset.*/
    dataset = H5Dcreate(root_id, "points", H5T_NATIVE_DOUBLE_g,
			dataspace, H5P_DEFAULT);
    if(dataset >= 0) {
      fld_points = (double *)malloc(2*fld_npoints*sizeof(double));
      for (ii=0; ii<fld_dims[0]; ii++) {
	fld_points[2*ii] = xcoord[ii];
	fld_points[2*ii+1] = ycoord[0];
      }
      for (ii=0; ii<fld_dims[1]; ii++) {
	fld_points[2*(ii+fld_dims[0])] = xcoord[0];
	fld_points[2*(ii+fld_dims[0])+1] = ycoord[ii];
      }
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace, H5P_DEFAULT,
	       fld_points);
      H5Dclose(dataset);	/* close dataset */
    }
    H5Sclose(dataspace);
  }

  /* Close the root group */
  H5Gclose(root_id);

  /* Restore standard HDF5 error function */
  H5Eset_auto(olderrfunc, stderr);

  /* Close the file */
  status = H5Fclose(file_id);

  return 0;
}

/************************************************************/
int writehdf5_(int *islice, int *ivar, double *data, int *ny, int *nx,
	       char *vname, long int namelen)
{
  char gname[64];
  int  nvar;
  static int nsteps=0;
  herr_t     status;
  hid_t      att_id, dataset, dataspace, file_id, timegroup, vargroup;
  hsize_t    dims[2];
  H5E_auto_t olderrfunc;

  /* Turn off HDF5 error messages */
  H5Eget_auto(&olderrfunc, NULL);
  H5Eset_auto(NULL, NULL);

  /* Open the HDF5 file for reading and writing */
  if ((file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT)) < 0) {
    fprintf(stderr, "Failed to open HDF5 file %s.\n", filename);
    return (int)file_id;
  }

  /* Determine the group name for this time slice */
  sprintf(gname, "time_node_data[%d]", *islice);
  /*sprintf(gname, "/");*/

  /* Try to open the group as if it existed */
  if ((timegroup = H5Gopen(file_id, gname)) < 0) {
    /* If not, create the group */
    if ((timegroup = H5Gcreate(file_id, gname, 0)) < 0) {
      fprintf(stderr, "Failed to create time group %s.\n", gname);
      H5Fclose(file_id);
      return (int)timegroup;
    }
    
    /* Initialize nnode_date attribute */
    write_int_attr(timegroup, "nnode_data", 0);

    /* Update number of time slices */
    nsteps++;
  }

  /* Update number of variables in time slice */
  att_id = H5Aopen_name(timegroup, "nnode_data");
  status = H5Aread(att_id, H5T_NATIVE_INT_g, &nvar);
  nvar++;
  status = H5Awrite(att_id, H5T_NATIVE_INT_g, &nvar);
  H5Aclose(att_id);

  /* Determine the group name for this variable */
  sprintf(gname, "node_data[%d]", *ivar);

  /* Create the variable group */
  if ((vargroup = H5Gcreate(timegroup, gname, 0)) < 0) {
    fprintf(stderr, "Failed to create group %s.\n", gname);
    H5Gclose(timegroup); H5Fclose(file_id);
    return (int)vargroup;
  }

  /* Vector length (=1 for scalars) */
  write_int_attr(vargroup, "veclen", 1);

  /* Save the variable name as an attribute */
  vname[namelen] = 0;
  write_string_attr(vargroup, "labels", vname);

  /* Units go here */
  /* write_string_attr(vargroup, "units", "furlongs per fortnight"); */

  /* Create a two-dimensional data space for the variable data */
  dims[0] = *nx;  dims[1] = *ny;
  if ((dataspace = H5Screate_simple(2, dims, NULL)) < 0) {
    fputs("Failed to create 2D variable dataspace.\n", stderr);
    return (int)dataspace;
  }

  /* Create the data set to hold the data */
  if ((dataset = H5Dcreate(vargroup, "values", H5T_NATIVE_DOUBLE, dataspace,
			   H5P_DEFAULT)) < 0) {
    fputs("Failed to create variable dataset.\n", stderr);
    H5Sclose(dataspace); H5Gclose(vargroup);
    H5Gclose(timegroup); H5Fclose(file_id);
    return (int)dataset;
  }

  /* Write the data values to the dataset in the file */
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, dataspace,
		    H5P_DEFAULT, (void *)data);

  /* Close the dataset and data space */
  status = H5Dclose(dataset);
  status = H5Sclose(dataspace);

  /* Close the variable and time slice groups and the data file */
  status = H5Gclose(vargroup);
  status = H5Gclose(timegroup);

  /* Update the nsteps info */
  timegroup = H5Gopen(file_id, "/");
  write_int_attr(timegroup, "nsteps", nsteps);
  status = H5Gclose(timegroup);

  /* Close the file */
  status = H5Fclose(file_id);

  /* Restore standard HDF5 error function */
  H5Eset_auto(olderrfunc, stderr);

  return 0;
}

/************************************************************/
int writehdf5scalar_(int *islice, int *ivar, double *value,
		    char *vname, long int namelen)
{
  char gname[64];
  herr_t     status;
  hid_t      file_id, timegroup, vargroup;
  H5E_auto_t olderrfunc;

  /* Turn off HDF5 error messages */
  H5Eget_auto(&olderrfunc, NULL);
  H5Eset_auto(NULL, NULL);

  /* Open the HDF5 file for reading and writing */
  if ((file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT)) < 0) {
    fprintf(stderr, "Failed to open HDF5 file %s.\n", filename);
    return (int)file_id;
  }

  /* Determine the group name for this time slice */
  sprintf(gname, "time_node_data[%d]", *islice);

  /* Open the group */
  if ((timegroup = H5Gopen(file_id, gname)) < 0) {
    fprintf(stderr, "Failed to open time group %s.\n", gname);
    H5Fclose(file_id);
    return (int)timegroup;
  }

  /* Determine the group name for this variable */
  sprintf(gname, "node_data[%d]", *ivar);

  /* Open the group */
  if ((vargroup = H5Gopen(timegroup, gname)) < 0) {
    fprintf(stderr, "Failed to open variable group %s.\n", gname);
    H5Gclose(timegroup); H5Fclose(file_id);
    return (int)vargroup;
  }

  /* Write the information to an attribute */
  vname[namelen] = 0;
  write_double_attr(vargroup, vname, *value);
  
  /* Close the variable and time slice groups and the data file */
  status = H5Gclose(vargroup);
  status = H5Gclose(timegroup);
  status = H5Fclose(file_id);

  /* Restore standard HDF5 error function */
  H5Eset_auto(olderrfunc, stderr);

  return 0;
}

/************************************************************/
int writehdf5time_(int *islice, double *time)
{
  char gname[64];
  hid_t file_id, timegroup;
  herr_t     status;
  H5E_auto_t olderrfunc;

  /* Open the HDF5 file for reading and writing */
  if ((file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT)) < 0) {
    fprintf(stderr, "Failed to open HDF5 file %s.\n", filename);
    return (int)file_id;
  }

  /* Turn off HDF5 error messages */
  H5Eget_auto(&olderrfunc, NULL);
  H5Eset_auto(NULL, NULL);

  /* Determine the group name for this time slice */
  sprintf(gname, "time_node_data[%d]", *islice);

  /* Try to open the group as if it existed */
  if ((timegroup = H5Gopen(file_id, gname)) < 0) {
    /* If not, create the group */
    if ((timegroup = H5Gcreate(file_id, gname, 0)) < 0) {
      fprintf(stderr, "Failed to create time group %s.\n", gname);
      H5Fclose(file_id);
      return (int)timegroup;
    }
  }

  /* Write the time value as an attribute */
  write_double_attr(timegroup, "time", *time);

  /* Close the group and file */
  status = H5Gclose(timegroup);
  status = H5Fclose(file_id);

  /* Restore standard HDF5 error function */
  H5Eset_auto(olderrfunc, stderr);

  return 0;
}

void __ctype_b(){}
