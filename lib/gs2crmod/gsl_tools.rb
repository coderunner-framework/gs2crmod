###################################
# GS2 CodeRunner Module: GSL Tools
# 
# Various useful methods for manipulating gsl data.
#
###################################

module GSL
	def self.cache
		@cache ||= {}
		@cache
	end
end


class GSL::Vector
	def pieces(no_pieces)
		ans = []
		piece_sizes  = []

		for i in 0...no_pieces 
			ans.push []; piece_sizes[i] = 0
		end
		for j in 0...size 
			piece_sizes[j % no_pieces] += 1
		end
# 		p ans, piece_sizes
		accum = 0
		piece_sizes.each_with_index do |piece_size, piece|
			ans[piece] = self.subvector(accum, piece_size)
			accum += piece_size
		end
		return ans
	end
	
	def from_box_order
		size = self.size
		if size < 3
			return self.dup
		end
		v1, v2 = self.subvector(0, (size+2)/2), self.subvector((size+2)/2, (size-1)/2)
		return v2.connect(v1)
	end
	def to_box_order
		size = self.size
		if size < 3
			return self.dup
		end
		v1, v2 = self.subvector(0, (size-1)/2), self.subvector((size-1)/2, (size+2)/2)
		return v2.connect(v1)
	end

end

       

class GSL::Matrix::Complex

	def self.re_im(re, im)
		raise "Shape of real and imaginary matrices must match" unless re.shape == im.shape
		rows, cols = re.shape
		mat = alloc(rows, cols)
		for i in 0...rows
			for j in 0...cols
				mat[i,j] = GSL::Complex.alloc([re[i,j], im[i,j]])
			end
		end
		return mat
	end
	
	
	def backward_cols_c2c(normalise = false)
		gm = self.dup
		rows, cols = gm.shape
		table = GSL.cache[[:fft_table, :complex, rows]] ||= GSL::FFT::ComplexWavetable.alloc(rows)
		work = GSL.cache[[:fft_work, :complex, rows]] ||= GSL::FFT::ComplexWorkspace.alloc(rows)
		for i in 0...cols
			vec = gm.col(i)
			vec.backward!(table, work)
			for j in 0...rows
				gm[j,i] =  vec[j]
			end
		end
		gm = gm / rows if normalise
		gm
	end

	def backward_rows_c2c(normalise = false)
		gm = self.dup
		rows, cols = gm.shape
		table = GSL.cache[[:fft_table, :complex, cols]] ||= GSL::FFT::ComplexWavetable.alloc(cols)
		work = GSL.cache[[:fft_work, :complex, cols]] ||= GSL::FFT::ComplexWorkspace.alloc(cols)
		for i in 0...rows
			vec = gm.row(i)
			vec.backward!(table, work)
			for j in 0...cols
				gm[i,j] =  vec[j]
			end
		end
		gm = gm / cols if normalise
		gm
	end

	def forward_cols_c2c
		gm = self.dup
		rows, cols = gm.shape
		table = GSL.cache[[:fft_table, :complex, rows]] ||= GSL::FFT::ComplexWavetable.alloc(rows)
		work = GSL.cache[[:fft_work, :complex, rows]] ||= GSL::FFT::ComplexWorkspace.alloc(rows)
		for i in 0...cols
			vec = gm.col(i)
			vec.forward!(table, work)
			for j in 0...rows
				gm[j,i] =  vec[j]
			end
		end
		gm
	end

	def backward_rows_cc2r(normalise = false)
		gm = self.dup
		rows, cols = gm.shape
	# 	if cols%2 == 0
	# 		newcols = cols*2
	# 	else
	# 		newcols = cols*2 - 1		
	# 	end
		was_even = rows.times.inject(true) do |bool, i|
			bool and (gm[i, cols - 1].imag == 0.0)
		end
	# 	ep was_even
			
		if was_even
			newcols = cols * 2 - 2
		else
			newcols = cols * 2 - 1
		end
		gm_re = GSL::Matrix.alloc(rows, newcols)
			
		table = GSL.cache[[:fft_table, :real, newcols]] ||= GSL::FFT::RealWavetable.alloc(newcols)
		work = GSL.cache[[:fft_work, :real, newcols]] ||= GSL::FFT::RealWorkspace.alloc(newcols)
		row = GSL::Vector::Complex.alloc(cols)
	# 	p rows
		for i in 0...rows
	# 		p i
	# 		row = gm.row(i)
			(0...cols).each{|j| row[j] = gm[i,j]}
			if was_even
				vec = row.concat(row.subvector(1, row.size - 2).reverse.conjugate) if cols > 2
			else
				vec = row.concat(row.subvector(1, row.size - 1).reverse.conjugate) if cols > 1
			end
			vec.backward!(table, work)
			for j in 0...newcols
				gm_re[i,j] =  vec[j].real
			end
		end
		gm_re = gm_re / newcols.to_f if normalise
		gm_re
	end

end

class GSL::Matrix

	def move_rows_from_box_order
		rows, cols = self.shape
		gm1, gm2 = self.view(0,0, (rows + 1)/2, cols), self.view((rows + 1)/2, 0,  (rows - 1)/2, cols)
		return gm2.vertcat(gm1)
	end
	def move_cols_from_box_order
		rows, cols = self.shape
		gm1, gm2 = self.view(0,0, rows, (cols + 1)/2), self.view(0, (cols + 1)/2, rows, (cols - 1)/2)
		return gm2.horzcat(gm1)
	end
	
	
def forward_rows_r2cc
	gm = self.dup
	rows, cols = gm.shape
	if cols%2 == 0
		newcols = (cols+2)/2  
	else
		newcols = (cols+1)/2		
	end
	gm_cc = GSL::Matrix::Complex.alloc(rows, newcols)
		
	table = GSL.cache[[:fft_table, :real, cols]] ||= GSL::FFT::RealWavetable.alloc(cols)
	work = GSL.cache[[:fft_work, :real, cols]] ||= GSL::FFT::RealWorkspace.alloc(cols)
	row = GSL::Vector.alloc(cols)
	for i in 0...rows
		(0...cols).each{|j| row[j] = gm[i,j]}
# 		p i
# 		row = gm.get_row(i)
# 		vec_out = GSL::Vector::Complex.alloc(row.size)
# 		(0...row.size).each{|j| vec[j] = row[j]}
# 		if cols%2 == 0
# 			vec_out = vec.concat(vec.subvector(1, vec.size - 2).reverse.conjugate)
# 		else
# 			vec_out = vec.concat(vec.subvector(1, vec.size - 1).reverse.conjugate)
# 		end
		view = row.forward(table, work).halfcomplex_to_complex
		for j in 0...newcols
			gm_cc[i,j] =  view[j]
		end
	end
	return gm_cc
end
	
end


class NArray
  def expand(*new_shape, empty_value)
    na = NArray.new(self.typecode,*new_shape)
		na[true] = empty_value 
    range = self.shape.map{|n| 0...n}
    na[*range] = self
    return na
  end
end
