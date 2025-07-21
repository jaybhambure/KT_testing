
	void get_cells_p(int L_coord, int *cellp)
	{	
		
		*cellp = (L_coord + 1)%L;

	}


	void get_cells_m(int L_coord, int *cellm)
	{	
		
		*cellm = (L_coord - 1 + L)%L;

	}




