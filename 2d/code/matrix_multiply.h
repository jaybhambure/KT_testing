void multiply(double Lambda[dim][dim], double a[dim][dim], double *result)
{

double product1[dim][dim] = {};
double product2[dim][dim] = {};

for (int i = 0; i < dim; i++)
{
	for (int j = 0; j < dim; j++)
	{
		for (int k = 0; k < dim; k++)
		{
			product1[i][j] += Lambda[i][k] * a[k][j];
		}
	}
}	


for (int i = 0; i < dim; i++)
{
	for (int j = 0; j < dim; j++)
	{
		for (int k = 0; k < dim; k++)
		{
			product2[i][j] += product1[i][k] * Lambda[j][k];
		}
	}
}
	
/*
for (int i = 0; i < dim; i++)
{
	for (int j = 0; j < dim; j++)
	{
		cout<<Lambda[i][j]<<" ";
	}
	cout<<endl;
}

cout<<endl<<endl;
*/


for (int i = 0; i < dim; i++)
{
	for (int j = 0; j < dim; j++)
	{
		int L_coord = i + dim*j;
		result[L_coord] = product2[i][j];
		//cout<<product2[i][j]<<" ";
	}
	//cout<<endl;
}

//cout<<endl<<endl;


}
