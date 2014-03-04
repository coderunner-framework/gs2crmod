#include "gs2crmod_ext.h"
#include <math.h>
#include <string.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline.h>

VALUE gs2crmod_tensor_complexes_field_gsl_tensor_complex_2(VALUE self, VALUE options)
{
	VALUE field, field_complex, field_narray;
	VALUE field_complex_narray;
	VALUE cgsl_tensor_complex;
	VALUE ccomplex;
	VALUE shape;
	int *c_shape;
	int i, j, k;

	ccomplex = RGET_CLASS_TOP("Complex");
	field = RFCALL_11("field_gsl_tensor", options);
	cgsl_tensor_complex = RGET_CLASS(cgsl, "TensorComplex");
	shape = rb_funcall(field, rb_intern("shape"), 0);
	/*rb_p(shape);*/

	CR_INT_ARY_R2C_STACK(shape, c_shape);
	field_complex = rb_funcall3(cgsl_tensor_complex, rb_intern("alloc"), 3, RARRAY_PTR(shape));
	field_narray = RFCALL_10_ON(field, "narray");
	field_complex_narray = RFCALL_10_ON(field_complex, "narray");
	/*rb_p(field_complex);*/
	for (i=0;i<c_shape[0];i++)
	  for (j=0;j<c_shape[1];j++)
			for (k=0;k<c_shape[2];k++)
			{

				/*printf("%d %d %d\n", i, j, k);*/
				VALUE cmplex, re, im;
				re = rb_funcall(
								field_narray,
								rb_intern("[]"),
								4,
								INT2FIX(0),
								INT2FIX(k), INT2FIX(j), INT2FIX(i)
							); 
				/*printf("%d %d %d re\n", i, j, k);*/
				im = rb_funcall(
								field_narray,
								rb_intern("[]"),
								4,
								INT2FIX(1),
								INT2FIX(k), INT2FIX(j), INT2FIX(i)
							);
				/*printf("%d %d %d im\n", i, j, k);*/
				cmplex = 
					 	rb_funcall(
							ccomplex,
							rb_intern("rectangular"),
							2, re, im
						);
				/*printf("%d %d %d again\n", i, j, k);*/
				rb_funcall(
						field_complex_narray, 
						rb_intern("[]="), 
						4, 
						INT2FIX(k), INT2FIX(j), INT2FIX(i),
						cmplex
					);

			}
	/*rb_p(field);*/



		/*field_complex = */
	return field_complex;
}

VALUE gs2crmod_tensor_field_gsl_tensor(VALUE self, VALUE options)
{
	VALUE cgsl_tensor;
	VALUE field, field_narray;
	VALUE field_real_space, field_real_space_narray;
	VALUE field_real_space_new, field_real_space_new_narray;
	VALUE shape, shape_real_space, shape_real_space_new;
	VALUE workspacex, workspacey;
	VALUE ccomplex, cgsl_complex;
	VALUE re, im, cmplex;
	VALUE *ary_ptr;
	VALUE range_x, range_y;
	int shape0;
	double* c_field, *c_field_real_space;
	int* c_shape, *c_shape_real_space;
	int * c_shape_real_space_new;
	int i,j,k,m;
	int kint,km,kold;
	double frac, interp_value;

	ccomplex = RGET_CLASS_TOP("Complex");
	cgsl_complex = RGET_CLASS(cgsl, "Complex");
	cgsl_tensor = RGET_CLASS(cgsl, "Tensor");
	/*cgsl_vector = RGET_CLASS(cgsl, "Vector");*/

	Check_Type(options, T_HASH);
	if(RTEST(CR_HKS(options, "field"))){
		field = CR_HKS(options, "field");
		CR_CHECK_CLASS(field, cgsl_tensor);
	}
	else {
	 field = RFCALL_11("field_gsl_tensor", options);
	}
	field_narray = RFCALL_10_ON(field, "narray");
	shape = RFCALL_10_ON(field, "shape");



	workspacex = RFCALL_11_ON(cgsl_vector_complex, "alloc", RARRAY_PTR(shape)[1]);
	shape0 = NUM2INT(RARRAY_PTR(shape)[0]);
	workspacey = RFCALL_11_ON(cgsl_vector, "alloc", INT2NUM(shape0 * 2 - 2 + shape0%2));

	field_real_space = rb_funcall(
			cgsl_tensor,
			rb_intern("alloc"),
			3,
			RFCALL_10_ON(workspacey, "size"),
			RARRAY_PTR(shape)[1],
			RARRAY_PTR(shape)[2]
			);
	field_real_space_narray = RFCALL_10_ON(field_real_space, "narray");
	shape_real_space = RFCALL_10_ON(field_real_space, "shape");
	
	CR_INT_ARY_R2C_STACK(shape, c_shape);
	CR_INT_ARY_R2C_STACK(RFCALL_10_ON(field_real_space, "shape"), c_shape_real_space);

	c_field = ALLOC_N(double, c_shape[0] * c_shape[1] * c_shape[2] * c_shape[3]);

	/*printf("Allocated stuff\n");*/
	for (j=0; j<c_shape[2]; j++) /*theta loop*/
	{
		for (i=0; i<c_shape[0]; i++) /*ky loop*/
		{
			for (k=0; k<c_shape[1]; k++) /*kx loop*/
			{
				/*First copy the field onto the 
				 * x workspace, Fourier transform
				 * the workspace and copy it onto
				 * the c array*/
				re = rb_funcall(
								field_narray,
								rb_intern("[]"),
								4,
								INT2FIX(0),
								INT2FIX(j), INT2FIX(k), INT2FIX(i)
							); 
				/*printf("%d %d %d re\n", i, j, k);*/
				im = rb_funcall(
								field_narray,
								rb_intern("[]"),
								4,
								INT2FIX(1),
								INT2FIX(j), INT2FIX(k), INT2FIX(i)
							);
				/*printf("%d %d %d im\n", i, j, k);*/
				cmplex = 
					 	rb_funcall(
							cgsl_complex,
							rb_intern("alloc"),
							2, re, im
						);
				/*printf("%d %d %d again\n", i, j, k);*/
				/*printf("Made complex\n");*/
				rb_funcall(
					workspacex,
					rb_intern("[]="),
					2,
					INT2FIX(k),
					cmplex);
				/*printf("Added complex\n");*/

			}
			/*printf("Made complex vector\n");*/
			workspacex = RFCALL_10_ON(workspacex, "backward");
			/*printf("Done x transform\n");*/
			for (k=0;k<c_shape[1];k++){
				cmplex = RFCALL_11_ON(workspacex,"[]",INT2FIX(k));
				c_field[
					j*c_shape[0]*c_shape[1]*2 +
					i*c_shape[1]*2 +
					k*2 
					] = NUM2DBL(RFCALL_10_ON(cmplex,"real"));
				c_field[
					j*c_shape[0]*c_shape[1]*2 +
					i*c_shape[1]*2 +
					k*2 + 1
					] = NUM2DBL(RFCALL_10_ON(cmplex,"imag"));

			}
		}
		/*printf("Done x\n");*/
		/* Now copy onto the y workspace,
		 * Fourier transform and copy onto 
		 * the second c array*/
		for (k=0;k<c_shape[1];k++)
		{
		  m=0;	
		  for (i=0;i<c_shape[0];i++)
			{
				rb_funcall(
						workspacey,
						rb_intern("[]="),
						2,
						INT2FIX(m),
						rb_float_new(c_field[
								j*c_shape[0]*c_shape[1]*2 +
								i*c_shape[1]*2 +
								k*2 
								])
						);
				m++;
				/* We are converting a complex array 
				 * to a half complex array in
				 * preparation for transforming 
				 * it to real space, and so there
				 * are one or two elements we leave
				 * unfilled.*/
				if (i==0 || (c_shape[0]%2==0 && i == c_shape[0]/2 + 1)) continue;
				rb_funcall(
						workspacey,
						rb_intern("[]="),
						2,
						INT2FIX(m),
						rb_float_new(c_field[
								j*c_shape[0]*c_shape[1]*2 +
								i*c_shape[1]*2 +
								k*2 + 1
								])
						);
				m++;
			}
			/*printf("Made y vector\n");*/
			/*rb_p(workspacey);*/
			/*rb_p(RFCALL_10_ON(workspacey, "class"));*/
				workspacey = RFCALL_10_ON(workspacey, "backward");
				/*printf("Done y transform\n");*/
			for (i=0;i<FIX2INT(RFCALL_10_ON(workspacey, "size"));i++)
			{
				rb_funcall(
						field_real_space_narray,
						rb_intern("[]="), 4,
						INT2FIX(j), INT2FIX(k), INT2FIX(i),
						RFCALL_11_ON(workspacey, "[]", INT2FIX(i))
						);
			}
		}
	}
	/*printf("HERE!\n");*/
  range_x = CR_RANGE_INC(
				(RTEST(CR_HKS(options, "ymin")) ?
					FIX2INT(CR_HKS(options, "ymin")) :
			 		0),
				(RTEST(CR_HKS(options, "ymax")) ?
					FIX2INT(CR_HKS(options, "ymax")) :
					c_shape_real_space[0]-1)
				);
	/*rb_p(range_x);*/
	range_y = CR_RANGE_INC(
				(RTEST(CR_HKS(options, "xmin")) ?
					FIX2INT(CR_HKS(options, "xmin")) :
			 		0),
				(RTEST(CR_HKS(options, "xmax")) ?
					FIX2INT(CR_HKS(options, "xmax")) :
					c_shape_real_space[1]-1)
				);
	/*printf("Made ranges\n");*/
	/*rb_p(rb_funcall(field_real_space, rb_intern("[]"), 3, INT2FIX(0), INT2FIX(2), INT2FIX(3)));*/
	field_real_space = rb_funcall(
			field_real_space,
			rb_intern("[]"), 3,
			range_x,
			range_y,
			Qtrue);
	/*rb_p(rb_funcall(field_real_space, rb_intern("[]"), 3, INT2FIX(0), INT2FIX(2), INT2FIX(3)));*/
	/*printf("SHORTENED!\n");*/
	if (RTEST(CR_HKS(options, "interpolate_theta")))
	{
		shape_real_space_new = RFCALL_10_ON(shape_real_space, "dup");
		kint = NUM2INT(CR_HKS(options, "interpolate_theta"));
		rb_funcall(
				shape_real_space_new,
				rb_intern("[]="), 2, INT2FIX(-1),
				INT2FIX((c_shape_real_space[2] - 1)*kint+1)
				);
		CR_INT_ARY_R2C_STACK(shape_real_space_new, c_shape_real_space_new); 
		ary_ptr = RARRAY_PTR(shape_real_space_new);
		field_real_space_new = rb_funcall(
				cgsl_tensor,
				rb_intern("float"),
				3,
				ary_ptr[0], ary_ptr[1], ary_ptr[2]);
		field_real_space_new_narray = RFCALL_10_ON(field_real_space_new, "narray");
		/*printf("Allocated");*/
		/*rb_p(shape_real_space_new);*/
		for (i=0;i<c_shape_real_space_new[0];i++)
		for (j=0;j<c_shape_real_space_new[1];j++)
		{
			/*printf("i %d j %d k %d\n", i, j, c_shape_real_space_new[2]-1); */
			rb_funcall(
					field_real_space_new_narray,
					rb_intern("[]="),
					4,INT2FIX(c_shape_real_space_new[2]-1),
					INT2FIX(j),INT2FIX(i),
					rb_funcall(
						field_real_space_narray,
						rb_intern("[]"), 3,
						INT2FIX(c_shape_real_space[2]-1),
						INT2FIX(j),INT2FIX(i)
						)
					);
			for (k=0;k<c_shape_real_space_new[2]-1;k++)
			{
				/*printf("i %d j %d k %d\n", i, j, k); */
				km = k%kint;
				frac = (double)km/(double)kint;
				kold = (k-km)/kint;
				interp_value = NUM2DBL(
						rb_funcall(
							field_real_space_narray,
							rb_intern("[]"),3,
							INT2FIX(kold),INT2FIX(j),INT2FIX(i)
							)
						)*(1.0-frac) + NUM2DBL(
						rb_funcall(
							field_real_space_narray,
							rb_intern("[]"),3,
							INT2FIX(kold+1),INT2FIX(j),INT2FIX(i)
							)
						) * frac;
				rb_funcall(
						field_real_space_new_narray,
						rb_intern("[]="), 4,
						INT2FIX(k),INT2FIX(j),INT2FIX(i),
						rb_float_new(interp_value)
						);
				/*if (i==0 && j==2 && k==3){*/
					/*printf("frac %f\n", frac);*/
				/*}*/
			}
		}

		field_real_space = field_real_space_new;
	}
	/*rb_p(shape_real_space_new);*/
	return field_real_space;
}

inline int arr_index(int ix, int iy, int ith, int it, int *c_shape, int t_size)
{
    return ix*c_shape[0]*c_shape[2]*t_size + iy*c_shape[0]*t_size + ith*t_size + it;
}

//Calculates normalized correlation function of 2 time series and stores result
void correlate_norm(double *A, double *B, double *C, int t_size)
{
    //printf("%lf, %lf, %d\n", A[0], B[0], size);

    int i, j, jmax;
    double temp_sum, norm_a=0., norm_b=0.;

    //printf("%lf, %lf, %lf, %d\n", A[0], B[0], C[0], t_size);
    
    //Sum over appropriate terms to get vector of correlations as function of time delay
    //Split into two loops instead of trying to be too smart with the indices.
    for(i=0; i<=t_size; i++){
        temp_sum=0;
        for(j=0; j<=i; j++)
        {
            temp_sum += A[j] * B[t_size - i - 1 + j];
        }
        C[i] = temp_sum;
    }
    //Second half
    for(i=t_size+1; i<2*t_size-1; i++){
        jmax = 2*t_size - 2 - i;
        temp_sum=0;
        for(j=0; j<=2*t_size-2-i; j++){
            temp_sum += A[t_size - jmax - 1 + j] * B[j];
        }
        C[i] = temp_sum;
    }

//    for(i=0; i<t_size; i++)
//       printf("%lf, ", A[i]);
//    printf("\n\n");
//    for(i=0; i<2*t_size-1; i++)
//        printf("%lf, ", C[i]);
//    printf("\n\ntemp_sum=%lf, Cmax=%lf\n", temp_sum, C[50]);

    //Normalize each array with its zero time delay correlation value
    for(i=0; i<t_size; i++){
        norm_a += A[i]*A[i];
        norm_b += B[i]*B[i];
    }
    //Final answer written to C array
    for(i=0; i<2*t_size-1; i++)
        C[i] = C[i]/sqrt(norm_a * norm_b);
}

VALUE gs2crmod_tensor_field_correlation_gsl_tensor(VALUE self, VALUE options)
{
    VALUE field_graphkit, shape, datakit, datakit_array;
    VALUE x_narray, y_narray, z_narray, field_narray;
    VALUE t_gsl_vector;
    //Write some as floats to avoid flooring to wrong (lower) number later
    float *lx_bin, *ly_bin, *lz_bin, *lt_bin;
    float lxmin, lxmax, lymin, lymax, lzmin, lzmax, ltmin, ltmax, lx, ly, lz, *lt;
    double *x, *y, *z, *t, *field, *coarray, *corr_norm, amin;
    int *c_shape, tot_size, ith, ix, iy, it, t_size, first=1;
    int index, i, i2, i3, i4, tot_bins, *count;
    long int i1;

    printf("Starting correlation analysis\n");

    /*Find time steps*/
    /*self.gsl_vector('t')*/
    //t_gsl_vector = RFCALL_11("gsl_vector", rb_str_new2("t"));
    t_gsl_vector = rb_funcall(self, rb_intern("gsl_vector"), 2, rb_str_new2("t"), options);
    t_size = NUM2INT(RFCALL_10_ON(t_gsl_vector, "size"));
    t = ALLOC_N(double, t_size);
    for(it=0; it<t_size; it++)
    {
        t[it] = NUM2DBL(CR_TELMT_R1(t_gsl_vector, it));
    }

    //Read from options
	  if(RTEST(CR_HKS(options, "amin")))
        amin = NUM2DBL(CR_HKS(options, "amin"));
    else
        amin = 1.0;

    VALUE nbins_array=CR_HKS(options, "nbins_array");
    int *nbins;
	  if(RTEST(nbins_array) && RTEST(rb_obj_is_kind_of(nbins_array, RGET_CLASS_TOP("Array")))){
        CR_INT_ARY_R2C_STACK(nbins_array, nbins);
    }
    else
        rb_raise(RGET_CLASS_TOP("TypeError"), "Please specify nbins_array as a 4D array");

    //correlation_type = options[:correlation_type]
    //Test which correlation type is to be calculated using ruby string comparison in a proc
    VALUE test_proc = rb_eval_string("Proc.new {|options| case options[:correlation_type]; when 'perp'; 0; when 'par'; 1; when 'time'; 2; when 'full'; 3 ;else; raise 'Please specify correlation_type as a string (perp/par/time/full)'; end}");
    int corr_type = FIX2INT(RFCALL_11_ON(test_proc, "call", options));
    
    for(it=0; it<t_size; it++)
    {   
        /*options[:t_index] = it+1*/
        /*options.send("[]=", :t_index, it+1)*/

        //rb_funcall(options, rb_intern("[]="), ID2SYM(rb_intern("t_index")), Qnil);
        CR_HKS_SET(options, "t_index", INT2FIX(it));

        field_graphkit = RFCALL_11("field_real_space_graphkit", options);

        /*shape = rb_funcall(field_graphkit, rb_intern("shape"), 0);*/
        datakit_array = RFCALL_10_ON(field_graphkit, "data");
        datakit = CR_ELEMENT_ACCESS(datakit_array, INT2FIX(0));
        
        /*Access data: datakit.x.data.narray*/ 
        x_narray = RFCALL_10_ON(RFCALL_10_ON(RFCALL_10_ON(datakit, "x"), "data"), "narray");
        y_narray = RFCALL_10_ON(RFCALL_10_ON(RFCALL_10_ON(datakit, "y"), "data"), "narray");
        z_narray = RFCALL_10_ON(RFCALL_10_ON(RFCALL_10_ON(datakit, "z"), "data"), "narray");
        field_narray = RFCALL_10_ON(RFCALL_10_ON(RFCALL_10_ON(datakit, "f"), "data"), "narray");

        shape = RFCALL_10_ON(x_narray, "shape");
        CR_INT_ARY_R2C_STACK(shape, c_shape);

        if(first)
        {
            tot_size = c_shape[0]*c_shape[1]*c_shape[2]*t_size; 
            x = ALLOC_N(double, tot_size);
            y = ALLOC_N(double, tot_size);
            z = ALLOC_N(double, tot_size);
            field = ALLOC_N(double, tot_size);
            first=0;
        }


        /*Copy NArrays to C arrays*/
        for(ith=0; ith<c_shape[0]; ith++)
            for(ix=0; ix<c_shape[1]; ix++)
                for(iy=0; iy<c_shape[2]; iy++)
                {
                    index = arr_index(ix,iy,ith,it,c_shape,t_size);
                    x[index] = NUM2DBL(CR_TELMT_R3(x_narray, ith, ix, iy));
                    y[index] = NUM2DBL(CR_TELMT_R3(y_narray, ith, ix, iy));
                    z[index] = NUM2DBL(CR_TELMT_R3(z_narray, ith, ix, iy));
                    field[index] = NUM2DBL(CR_TELMT_R3(field_narray, ith, ix, iy));
                }
    }

    

    /******************************
     * GSL Interpolation of field *
     *        in time             *
     *****************************/
    //Read t interp from options (default is 100)
    int nt_reg;
	  if(RTEST(CR_HKS(options, "nt_reg")))
        nt_reg = NUM2INT(CR_HKS(options, "nt_reg"));
    else
        nt_reg = 100;
    int idx, ti;
    double *x_reg, *y_reg, *z_reg, *t_reg, *field_reg, delta_t_reg, *y1;
    tot_size = c_shape[0]*c_shape[1]*c_shape[2]*nt_reg; 
    //printf("Start interpolation, %d, %d, %d, %d, %d\n", c_shape[0], c_shape[1], c_shape[2], nt_reg, tot_size);
    delta_t_reg = (t[t_size-1]-t[0])/(nt_reg-1);

    y1 = ALLOC_N(double, t_size);
    x_reg = ALLOC_N(double, tot_size);
    y_reg = ALLOC_N(double, tot_size);
    z_reg = ALLOC_N(double, tot_size);
    t_reg = ALLOC_N(double, nt_reg);
    field_reg = ALLOC_N(double, tot_size);

    for(i1=0; i1<c_shape[0]*c_shape[1]*c_shape[2]*t_size; i1+=t_size)
    {
        //printf("%d, %d\n", i1, c_shape[0]*c_shape[1]*c_shape[2]*t_size);
        for (i2 = 0; i2 < t_size; i2++)
        {
            y1[i2] = field[i1+i2];
            //printf ("%d, %lf %lf\n", i2, t[i2], y1[i2]);
        }

        gsl_interp_accel *acc 
              = gsl_interp_accel_alloc ();
        gsl_spline *spline 
              = gsl_spline_alloc (gsl_interp_cspline, t_size);

        gsl_spline_init (spline, t, y1, t_size);

        //printf("Interpolated:\n");
        for (ti = 0; ti < nt_reg; ti++)
        {
            t_reg[ti] = t[0]+ti*delta_t_reg;
            idx = floor((i1)/t_size)*nt_reg;
            x_reg[idx+ti] = x[i1];
            y_reg[idx+ti] = y[i1];
            z_reg[idx+ti] = z[i1];
            field_reg[idx+ti] = gsl_spline_eval (spline, t_reg[ti], acc);
        }
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);

    }
    
    //Can now free the original pointers
    free(x); x=0;
    free(y); y=0;
    free(z); z=0;
    free(t); t=0;

    /*printf("interpolated:\n");
    for(i1=0; i1<c_shape[1]; i1++){
        index = arr_index(i1,0,0,0,c_shape,nt_reg);
        printf("%lf, ", x_reg[index]);
    }*/

    /**********************************************************************
     * First need to know max and min lengths and times in each dimension *
     * in order to define bins. These need to be specified before the loop*
     * since the binning will be done at each step since there is too much*
     * info to store.                                                     *
     * ********************************************************************/
    double lx_pos_min=20;
    i3=0;
    lxmin=0; lymin=0; lzmin=0;
    for(i1=0; i1<tot_size; i1+=nt_reg){
        //printf("%lf, %lf, %lf\n", x_reg[i1], y_reg[i1], z_reg[i1]);
        for(i2=0; i2<tot_size; i2+=nt_reg)
        {
            lx = x_reg[i2] - x_reg[i1];
            ly = y_reg[i2] - y_reg[i1];
            lz = z_reg[i2] - z_reg[i1];
            
            if(lx < lxmin)
                lxmin = lx;
            if(ly < lymin)
                lymin = ly;
            if(lz < lzmin)
                lzmin = lz;
        }
    }
    lxmax = -lxmin;
    lymax = -lymin;
    lzmax = -lzmin;
    ltmin = t_reg[0] - t_reg[nt_reg-1];
    ltmax = -ltmin;
    //printf("\nMin/Max = %lf, %lf\n", lxmin, lxmax);
    //printf("Min/Max = %lf, %lf\n", lymin, lymax);
    //printf("Min/Max = %lf, %lf\n", lzmin, lzmax);

    //printf("Bins\n");
    //Initialize the bin arrays
    tot_bins = nbins[0]*nbins[1]*nbins[2]*nbins[3];
    lx_bin = ALLOC_N(double, nbins[0]);
    ly_bin = ALLOC_N(double, nbins[1]);
    lz_bin = ALLOC_N(double, nbins[2]);
    lt_bin = ALLOC_N(double, nbins[3]);
    
    for(i1=0; i1<nbins[0]; i1++){
        lx_bin[i1] = lxmin + i1*(lxmax-lxmin)/(nbins[0]-1);
        //printf("lx_bin = %f\n ", lx_bin[i1]);
    }
    for(i1=0; i1<nbins[1]; i1++){
        ly_bin[i1] = lymin + i1*(lymax-lymin)/(nbins[1]-1);
        //printf("ly_bin = %f\n", ly_bin[i1]);
    }
    for(i1=0; i1<nbins[2]; i1++){
        lz_bin[i1] = lzmin + i1*(lzmax-lzmin)/(nbins[2]-1);
        //printf("lz_bin = %f\n", lz_bin[i1]);
    }
    for(i1=0; i1<nbins[3]; i1++)
        lt_bin[i1] = ltmin + i1*(ltmax-ltmin)/(nbins[3]-1);

    //Now define and initialize the coarray and count arrays
    coarray = ALLOC_N(double, tot_bins);
    count = ALLOC_N(int, tot_bins);
    corr_norm = ALLOC_N(double, (2*nt_reg-1));   //store correlate result before adding to coarray
    lt = ALLOC_N(double, (2*nt_reg-1));   //Can calculate lt before loops

    for(i1=0; i1<tot_bins; i1++){
        coarray[i1] = 0.;
        count[i1] = 0;
    }
    /******************************************************
     * Start looping and calculating correlation function *
     ******************************************************/
    //Can predefine the time separations
    for(i1=0; i1<2*nt_reg-1; i1++){
        if(i1<nt_reg)
            lt[i1] = t_reg[i1] - t_reg[nt_reg-1];
        else if(i1>nt_reg-1)
            lt[i1] = t_reg[i1-nt_reg] - t_reg[0] + delta_t_reg;
    }

    /* Now test which correlation function is to be calculated since full 4D correlation
     * is usually intractable. The type was read in at begininning and corr_type corresponds to:
     *
     * 0 : perpendicular only (lz = 0)
     * 1 : parallel only (lx = ly = 0)
     * 2 : time only (lx = lz = 0)
     * 3 : full correlation (may take very long) (all separations != 0)
     */
    float eps1 = 1e-5, eps2 = 1e5;  //define very small and very large numbers to test against
    float lx_test, ly_test, lz_test;
    switch (corr_type){
      case 0: //perp
        lx_test = eps2;
        ly_test = eps2;
        lz_test = eps1;
        break;
      case 1: //par
        lx_test = eps1;
        ly_test = eps1;
        lz_test = eps2;
        break;
      case 2: //time
        lx_test = eps1;
        ly_test = eps2;
        lz_test = eps1;
        break;
      case 3: //full
        lx_test = eps2;
        ly_test = eps2;
        lz_test = eps2;
        break;
    }

    //Start main loop for correlation function calculation
    for(i1=0; i1<tot_size; i1+=nt_reg){
        for(i2=i1; i2<tot_size; i2+=nt_reg)
        {
            //Calculate spatial and temporal separation
            lx = x_reg[i2] - x_reg[i1]; 
            ly = y_reg[i2] - y_reg[i1]; 
            lz = z_reg[i2] - z_reg[i1];

            if(lx<lx_test && ly<ly_test && lz<lz_test){
                //Calculate correlation function:
                //corr = correlate(field[i1], field[i2], answer, no of t pts)
                correlate_norm(&field_reg[i1], &field_reg[i2], &corr_norm[0], nt_reg);

                //Calculate appropriate bin (subtracting min ensures +ve idx)
                ix = floor((lx - lxmin) / (lx_bin[1] - lx_bin[0]));
                iy = floor((ly - lymin) / (ly_bin[1] - ly_bin[0]));
                ith = floor((lz - lzmin) / (lz_bin[1] - lz_bin[0]));
                //printf("(%lf, %lf, %lf, %lf)\n", lz, lzmin, lz_bin[1], lz_bin[0]);
                //printf("indices: (%d, %d, %d)\n", ix, iy, ith);

                //Loop over time calculate time bin and put into coarray
                for(i3=0; i3<2*nt_reg-1; i3++){
                    it = floor((lt[i3] - ltmin) / (lt_bin[1] - lt_bin[0]));
                    coarray[ix*nbins[1]*nbins[2]*nbins[3] + iy*nbins[2]*nbins[3] + ith*nbins[3] + it] += corr_norm[i3];
                    count[ix*nbins[1]*nbins[2]*nbins[3] + iy*nbins[2]*nbins[3] + ith*nbins[3] + it] += 1;
                }

                /**************************
                 * Repeat for negative    *
                 * lx, ly, lz, to avoid   *
                 * recalculating corr_norm*
                 **************************/

                lx = -lx; ly = -ly; lz = -lz;

                //Calculate appropriate bin (subtracting min ensures +ve idx)
                ix = floor((lx - lxmin) / (lx_bin[1] - lx_bin[0]));
                iy = floor((ly - lymin) / (ly_bin[1] - ly_bin[0]));
                ith = floor((lz - lzmin) / (lz_bin[1] - lz_bin[0]));
                //printf("- = {%d, %d, %d}\n", ix, iy, ith);

                //Loop over time calculate time bin and put into coarray
                for(i3=0; i3<2*nt_reg-1; i3++){
                    it = floor((-lt[i3] - ltmin) / (lt_bin[1] - lt_bin[0]));
                    coarray[ix*nbins[1]*nbins[2]*nbins[3] + iy*nbins[2]*nbins[3] + ith*nbins[3] + it] += corr_norm[i3];
                    count[ix*nbins[1]*nbins[2]*nbins[3] + iy*nbins[2]*nbins[3] + ith*nbins[3] + it] += 1;
                }
             }
         }
    } 
    //End main loop

    //Finally have to normalize coarray with count when count != 0
    //printf("Normalization: \n");
    for(i1=0; i1<tot_bins; i1++){
       if(count[i1]>0){
           coarray[i1] = coarray[i1]/count[i1];
       }
    }
    //printf("Finish normalization: \n");

    //Retun output to CR
	  VALUE cgsl_tensor, output_tensor;
    cgsl_tensor = RGET_CLASS(cgsl, "Tensor");
    //output_tensor = GSL::Tensor.alloc(nbins, ...)
    output_tensor = rb_funcall(cgsl_tensor, rb_intern("alloc"), 4, INT2FIX(nbins[0]), INT2FIX(nbins[1]), INT2FIX(nbins[2]), INT2FIX(nbins[3]));

    for(i1=0; i1<nbins[0]; i1++)
        for(i2=0; i2<nbins[1]; i2++)
            for(i3=0; i3<nbins[2]; i3++)
                for(i4=0; i4<nbins[3]; i4++){
                    CR_TELMT_R4_SET(output_tensor, i1, i2, i3, i4, rb_float_new(coarray[i1*nbins[1]*nbins[2]*nbins[3] + i2*nbins[2]*nbins[3] + i3*nbins[3] + i4]));
                }

    printf("Finished correlation analysis\n");
    return output_tensor; 
}

void Init_gs2crmod_ext()
{
	/*printf("HERE!!!");*/
	ccode_runner =  RGET_CLASS_TOP("CodeRunner");
	ccode_runner_gs2 =  rb_define_class_under(ccode_runner, "Gs2",
			RGET_CLASS(
				RGET_CLASS(ccode_runner, "Run"), 
				"FortranNamelist"
				)
			);

	ccode_runner_gs2_gsl_tensor_complexes = rb_define_module_under(ccode_runner_gs2, "GSLComplexTensors");
	rb_include_module(ccode_runner_gs2, ccode_runner_gs2_gsl_tensor_complexes);

	ccode_runner_gs2_gsl_tensors = rb_define_module_under(ccode_runner_gs2, "GSLTensors"); 
	rb_include_module(ccode_runner_gs2, ccode_runner_gs2_gsl_tensors);

	cgsl = RGET_CLASS_TOP("GSL");
	cgsl_vector = RGET_CLASS(cgsl, "Vector");
	cgsl_vector_complex = RGET_CLASS(cgsl_vector, "Complex");

	rb_define_method(ccode_runner_gs2_gsl_tensor_complexes, "field_gsl_tensor_complex_2", gs2crmod_tensor_complexes_field_gsl_tensor_complex_2, 1);
	rb_define_method(ccode_runner_gs2_gsl_tensors, "field_real_space_gsl_tensor", gs2crmod_tensor_field_gsl_tensor, 1);
	rb_define_method(ccode_runner_gs2_gsl_tensors, "field_correlation_gsl_tensor", gs2crmod_tensor_field_correlation_gsl_tensor, 1);

	/*rb_define_method(ccode_runner_ext, "hello_world", code_runner_ext_hello_world, 0);*/
}

