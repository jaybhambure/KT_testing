
	void get_cells_p(int L_coord, int *cellp, int e_x, int e_y)
	{	
		
		int j = L_coord/L;
		int i = L_coord - j*L;	

		int ip, jp;

		ip = (1 - e_x)*i + e_x * ( (i + 1)%L );

		jp = (1 - e_y)*j + e_y * ( (j + 1)%L );

		*cellp = ip + L*jp;			
	}


	void get_cells_m(int L_coord, int *cellm, int e_x, int e_y)
	{	
		
		int j = L_coord/L;
		int i = L_coord - j*L;	

		int im, jm;

		im = (1 - e_x)*i + e_x * ( (i - 1 + L)%L );

		jm = (1 - e_y)*j + e_y * ( (j - 1 + L)%L );

		*cellm = im + L*jm;	

	}




