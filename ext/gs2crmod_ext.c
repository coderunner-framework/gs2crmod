#include "gs2crmod_ext.h"

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

	/*rb_define_method(ccode_runner_ext, "hello_world", code_runner_ext_hello_world, 0);*/
}

