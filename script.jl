function systemMatrix(k1,k2,n)
        A = BandedMatrix(Eye(n),(2,2));
	A[band(-2)]=k1*ones(length(A[band(-2)])); # fill in the values for sub-diagonal 2
	A[band(-1)]=-8*k1*ones(length(A[band(-1)])); # fill in the values for sub-diagonal 1

	A[band(1)]=8*k1*ones(length(A[band(1)])); # fill in the values for super-diagonal 1
	A[band(2)]=-k1*ones(length(A[band(2)])); # fill in the values for super-diagonal 2

	AA = Matrix(A); # Convert to regular matrix type instead of banded matrix

	# For the first row, Dirichlet boundary condition
	AA[1,2]=0; AA[1,3]=0;
	
	# For the second and second-to-last rows, implement second-order central difference
	AA[2,1]= -k2; AA[2,3] = k2; AA[2,4] = 0;
	AA[end-1,end] = k2; AA[end-1,end-2]=-k2; AA[end-1,end-3] = 0;
	
	# For the last row, implement von Neumann boundary condition
	AA[end,end-1]=-1; AA[end,end-2] = 0;
	
	A= BandedMatrix(AA,(2,2)) # Convert back to banded matrix
	return A;
end