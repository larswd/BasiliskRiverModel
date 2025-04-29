/*
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on ASCII format. If one writes
one *.vts file per PID process this function may be combined with
output_pvtu_ascii() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
*/
double dt_start = 0.;
// Define alternative traversal order to match with the x ordering of vts (basilisk uses z ordering by default)

/*
Writes a structured vtk file (.vts) of a single layer.
*/
void output_vts_ascii_single_layer(FILE* fp, scalar* list, int layerid, bool displace_z, bool outputbreaking )
{

    
    int nthreads_ = omp_get_num_threads();
    omp_set_num_threads(1); 
    
        // MULTIGRID
    
 fputs("<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N, 0, 0);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, N, 0, 1);
#endif

#if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, 0);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, 1);
#endif

    // Loop over velocity data and store kinematics in cell vector stucture
    int layertemp = _layer;
    _layer = layerid;
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");
    foreach() {
#if dimension == 1
        fprintf(fp, "%g %g 0.\n", u.x[], w[]);
#endif
#if dimension == 2
        fprintf(fp, "%g %g %g\n", u.x[], u.y[], w[]);
#endif
    }
    fputs("\t\t\t\t </DataArray>\n", fp);


    // loop over all scalars in scalarlist and store values as cell data
    for (scalar s in list) {
        if (strcmp(s.name, "eta") == 0) {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            foreach() {
                if (h[] > dry) {
                    fprintf(fp, "%g\n", val(s));
                }
                else {
                    fprintf(fp, "nan\n");
                }
            }
            fputs("\t\t\t\t </DataArray>\n", fp);
        }
        else {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            foreach() {
                fprintf(fp, "%g\n", val(s));
            }
            fputs("\t\t\t\t </DataArray>\n", fp);
        }
        
    }

   

    fputs("\t\t\t </CellData>\n", fp);


    // Coordinates 
    fputs("\t\t\t <Points>\n", fp);
#if dimension == 1
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

    if (displace_z) {
        scalar z_trans = list[0];
        // set z coordinate to the value of the first scalar field in the provided list.
        double hh;
        foreach_vertex(serial) {
            hh = (z_trans[] + z_trans[-1]) / 2.;
            fprintf(fp, "%12.4f %12.4f 0.\n", x, hh);

        }
    }
    else {
        // set z coordinate to 0.
        foreach_vertex() {
            fprintf(fp, "%12.4f 0. 0.\n", x);
        }
    }
#endif
#if dimension == 2
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    
    if (displace_z) {
        scalar z_trans = list[0];
        // first do the bottom coordinates stored in zb
        double hh;
        foreach_vertex(serial) {
            hh = (z_trans[] + z_trans[-1] + z_trans[0, -1] + z_trans[-1, -1]) / 4.;
            fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh);
        }
    }
    else {
        // first do the bottom coordinates stored in zb
        foreach_vertex() {
            fprintf(fp, "%12.4f %12.4f 0.0\n", x, y);
        }
    }
#endif

    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);

    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    fputs("\t </StructuredGrid>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);

    _layer = layertemp;
    omp_set_num_threads(nthreads_);

}


void output_vts_ascii_all_layers(FILE* fp, scalar* list, int N)
{

    
    int nthreads_ = omp_get_num_threads();
    omp_set_num_threads(1); 
    
        // MULTIGRID
    
 fputs("<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, N, 0, nl);
#endif

#if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
#endif

    // Loop over velocity data and store kinematics in cell vector stucture
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");
    for(int i = nl-1; i >= 0; i--){
    foreach() {
#if dimension == 1
        fprintf(fp, "%g %g 0.\n", u.x[0,0,i], w[0,0,i]);
#endif
#if dimension == 2
        if (h[]> dry)
          fprintf(fp, "%g %g %g\n", u.x[0,0,i], u.y[0,0,i], w[0,0,i]);
        else
          fprintf(fp, "0 0 0\n");
#endif
    }
    }
    // Dummy text to get correct amount of layers
    foreach() {
        fprintf(fp, "0 0 0.\n");
    }
    fputs("\t\t\t\t </DataArray>\n", fp);


    // loop over all scalars in scalarlist and store values as cell data
    for (scalar s in list) {
        if (strcmp(s.name, "eta") == 0) {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            for(int i = nl-1; i >= 0; i--){
                foreach() {
                if (h[] > dry) {
                    fprintf(fp, "%g\n", s[0,0,i]);
                }
                else {
                    fprintf(fp, "nan\n");
                }
            }
        }
        fputs("\t\t\t\t </DataArray>\n", fp);
        }
        else {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            for(int i = nl-1; i >= 0; i--){
            foreach() {
                fprintf(fp, "%g\n", s[0,0,i]);
            }
            }
            foreach(){
                fprintf(fp, "0\n");
            }
            fputs("\t\t\t\t </DataArray>\n", fp);
        }   
    }   

    fputs("\t\t\t </CellData>\n", fp);
    // Coordinates 
    fputs("\t\t\t <Points>\n", fp);
#if dimension == 1
  fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
  double* zcorr = (double *)malloc((N+1)*sizeof(double));
  for (int j = 0; j <= N; j++){
    zcorr[j] = 0;
  }
  int k;
  for(int i = nl-1; i >= 0; i--){
    k = 0;
    foreach_vertex(serial) {
        fprintf(fp, "%12.4f %12.4f 0.\n", x,zcorr[k]);
        if (h[0,0,k] < 1e-3){
          zcorr[k] = zcorr[k] - h[-1,-1,i];
        } else{ 
          zcorr[k] = zcorr[k] - h[0,0,i];
        }
        k++;
    }
  }
  k = 0;
  foreach_vertex(serial){
    fprintf(fp, "%12.4f %12.4f 0.\n", x, zcorr[k]);
    k++;
  }
  free(zcorr);
#endif
#if dimension == 2
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
        
        double* zcorr = (double*)malloc((N+1)*(N+1)*sizeof(double));
        for (int j = 0; j < (N+1)*(N+1); j++){
            zcorr[j] = 0;
        }
        int j;
        for(int i = nl-1; i >= 0; i--){
            j = 0;
            foreach_vertex(serial) {
                int xshift, yshift;
                if ((x > X0) && (x < X0 + L0)){
                    xshift = 0;
                } else if (x < X0) {
                    xshift = 1;
                } else {
                    xshift = -1;
                }
                
                if ((y > Y0) && (y < Y0 + L0)){
                    yshift = 0;
                } else if (y < Y0){
                    yshift = 1;
                } else {
                    yshift = -1;
                }

                fprintf(fp, "%12.4f %12.4f %12.4f\n", x,y, zb[xshift, yshift] + h[xshift, yshift]);
                
                if (h[] < dry){
                    zcorr[j] = zcorr[j] - h[-1,-1,i];
                } else {
                    zcorr[j] = zcorr[j] - h[0,0,i];
                }
                j++;
            }
        }
        j = 0;
        foreach_vertex(serial){
                fprintf(fp, "%12.4f %12.4f %12.4f\n", x,y, zb[]);
                j++;
        }
        free(zcorr);
#endif

    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);

    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    fputs("\t </StructuredGrid>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);

    omp_set_num_threads(nthreads_);

}
