/*
 This code is written by Andrea Gobbi  <gobbi.andrea@mail.com> 2013

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "BiRewire.h"
#include <math.h>
static inline void loadBar(size_t x, size_t n, int r, int w)
{

    if(n<100)return;
    // Only update r times.
    if ( x % (n/r) != 0 ) return;

    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;

    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );

    // Show the load bar.
    for ( x=0; x<c; x++)
        printf("=");

    for ( x=c; x<w; x++)
        printf(" ");

    // ANSI verbose codes to go back to the
    // previous line and clear it.
    // printf("]\n33[F33[J");
    printf("]\r"); // Move to the first column

    ////fflush(stdout);
}
size_t inline min(size_t a,size_t b)
{
	if(a<b)
		return(a);
		return(b);
}





double inline unif_rand()
{
	return( (double)rand()/((double)RAND_MAX+1.0));
}







double inline similarity(short **m,short **n,size_t ncol,size_t nrow,size_t e)
{

	size_t num=0;
	size_t i, j;
	for(i=0;i<nrow;++i)
		for(j=0;j<ncol;++j)
			num += m[i][j] * n[i][j];
	return((double)num/(2.0*e-num));

}




size_t analysis(short **incidence,size_t ncol,size_t nrow,double **scores,size_t step,size_t max_iter,size_t verbose)
{



	size_t i,j,kk,n,rand1,rand2;
	size_t dim=max_iter+1;
	size_t *from;
  size_t *to;
  short **matrix;

	size_t a,b,c,d,e=0;
	size_t index=1;
  do matrix=(short **)calloc(nrow,sizeof(short*)); while(matrix==NULL);
	for(i=0;i<nrow;++i)
	{
		do	matrix[i]= (short*)calloc(ncol,sizeof(short)); while(matrix[i]==NULL);
    for(j=0;j<ncol;++j)
			{
				matrix[i][j]=incidence[i][j];
				e+=incidence[i][j];
			}
  }
	//initialization of score vector overwriting the original
	do	*scores=(double*)calloc(dim,sizeof(double));	while(scores==NULL);
	for(i=0;i<dim;(*scores)[i++]=0.0);
	(*scores)[0]=1.0;
 	do	from=(size_t*)calloc(e,sizeof(size_t));	while(from==NULL);
	do	to=(size_t*)calloc(e,sizeof(size_t));	  while(to==NULL);
 	kk=0;
	for(i=0;i<nrow;++i)
		for(j=0;j<ncol;++j)
			if(matrix[i][j]==1)
				{
					from[kk]=i;
					to[kk]=j;
					kk++;
				}
	time_t  tin,tfin;
	tin = time (NULL);
	//GetRNGstate();
	for(n=0;n<max_iter;n++)
	{
		//random rewiring
    if(verbose==1)
  		loadBar( n,  max_iter, 100,  50);
		rand1=(size_t) (unif_rand()*e);
		do	rand2=(size_t) (unif_rand()*e);	while (rand1==rand2);
		a=from[rand1];
		c=from[rand2];
		b=to[rand1];
		d=to[rand2];
		//printf("%d %d %d %d %d %d %d %d %d\n",e,a,b,c,d,rand1,rand2,t,n);
		if(a!=c && d!=b && incidence[a][d]==0 && incidence[c][b]==0)
			{
				to[rand1]=d;
				to[rand2]=b;
				incidence[a][d]=incidence[c][b]=1;
				incidence[a][b]=incidence[c][d]=0;

			}

			if(n%step==0)
			(*scores)[index++]=similarity(incidence,matrix, ncol,nrow, e);
	}
	tfin = time (NULL);
	//PutRNGstate();
	if(verbose==1)
    printf("DONE in %d seconds \n",-(tin-tfin));

	return (index-1);
}



size_t analysis_ex(short **incidence,size_t ncol,size_t nrow,double **scores,size_t step,size_t max_iter,size_t verbose, size_t MAXITER)
{



	size_t i,j,kk,n,rand1,rand2,t=0;
	size_t dim=max_iter+1;
	size_t *from;
  size_t *to;
  short **matrix;

	size_t a,b,c,d,e=0;
	size_t index=1;
  do matrix=(short **)calloc(nrow,sizeof(short*)); while(matrix==NULL);
	for(i=0;i<nrow;++i)
	{
		do	matrix[i]= (short*)calloc(ncol,sizeof(short)); while(matrix[i]==NULL);
    for(j=0;j<ncol;++j)
			{
				matrix[i][j]=incidence[i][j];
				e+=incidence[i][j];
			}
  }
	//initialization of score vector overwriting the original
	do	*scores=(double*)calloc(dim,sizeof(double));	while(scores==NULL);
	for(i=0;i<dim;(*scores)[i++]=0.0);
	(*scores)[0]=1.0;
 	do	from=(size_t*)calloc(e,sizeof(size_t));	while(from==NULL);
	do	to=(size_t*)calloc(e,sizeof(size_t));	  while(to==NULL);
 	kk=0;
	for(i=0;i<nrow;++i)
		for(j=0;j<ncol;++j)
			if(matrix[i][j]==1)
				{
					from[kk]=i;
					to[kk]=j;
					kk++;
				}
	time_t  tin,tfin;
	tin = time (NULL);
	//GetRNGstate();
	for(n=0;n<max_iter;t++)
	{
		//random rewiring
    if(verbose==1)
  		loadBar( n,  max_iter, 100,  50);
		rand1=(size_t) (unif_rand()*e);
		do	rand2=(size_t) (unif_rand()*e);	while (rand1==rand2);
		a=from[rand1];
		c=from[rand2];
		b=to[rand1];
		d=to[rand2];
		//printf("%d %d %d %d %d %d %d %d %d\n",e,a,b,c,d,rand1,rand2,t,n);
		if(a!=c && d!=b && incidence[a][d]==0 && incidence[c][b]==0)
			{
				to[rand1]=d;
				to[rand2]=b;
				incidence[a][d]=incidence[c][b]=1;
				incidence[a][b]=incidence[c][d]=0;
				n++;
			if(n%step==0)
			(*scores)[index++]=similarity(incidence,matrix, ncol,nrow, e);
			}
	if(t>MAXITER)
		{
			tfin = time (NULL);
			//PutRNGstate();
			if(verbose==1)
  	  printf("DONE in %d seconds \n",-(tin-tfin));

			printf("Reached the maximum number admissible of iterations!\n");
			return (index);
		}

	}
	tfin = time (NULL);
	//PutRNGstate();
	if(verbose==1)
    printf("DONE in %d seconds \n",-(tin-tfin));

	return (index-1);
}

size_t rewire_bipartite(short **matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose)
{
	size_t i,j,kk,n,e=0,rand1,rand2;
	size_t *from;
	size_t a,b,c,d;
	size_t *to;

	for(i=0;i<nrow;++i)
		for(j=0;j<ncol;++j)
            e+=matrix[i][j];

	do	from=(size_t*)calloc(e,sizeof(size_t));	while(from==NULL);

	do	to=(size_t*)calloc(e,sizeof(size_t));	while(to==NULL);
	kk=0;
	for(i=0;i<nrow;++i)
		for(j=0;j<ncol;++j)
			if(matrix[i][j]==1)
				{
					from[kk]=i;
					to[kk]=j;
					kk++;
				}
	//GetRNGstate();
	time_t  tin,tfin;
	tin = time (NULL);
	for(n=0;n<max_iter;n++)
	{
		//random rewiring
		if(verbose==1)
    	loadBar( n,  max_iter, 100,  50);
		rand1=(size_t) (unif_rand()*e);
		do  rand2=(size_t) (unif_rand()*e);while (rand1==rand2);
		a=from[rand1];
		c=from[rand2];
		b=to[rand1];
		d=to[rand2];
		if(a!=c && d!=b && matrix[a][d]==0 && matrix[c][b]==0)
			{
				to[rand1]=d;
				to[rand2]=b;
				matrix[a][d]=matrix[c][b]=1;
				matrix[a][b]=matrix[c][d]=0;

			}

	}
	//PutRNGstate();
	tfin = time (NULL);
	if(verbose==1)
    printf("DONE in %d seconds \n",-(tin-tfin));


return 0;

}


size_t rewire_bipartite_ex(short **matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose,size_t MAXITER)
{
	size_t i,j,kk,n,e=0,rand1,rand2,t=0;
	size_t *from;
	size_t a,b,c,d;
	size_t *to;

	for(i=0;i<nrow;++i)
		for(j=0;j<ncol;++j)
            e+=matrix[i][j];

	do	from=(size_t*)calloc(e,sizeof(size_t));	while(from==NULL);

	do	to=(size_t*)calloc(e,sizeof(size_t));	while(to==NULL);
	kk=0;
	for(i=0;i<nrow;++i)
		for(j=0;j<ncol;++j)
			if(matrix[i][j]==1)
				{
					from[kk]=i;
					to[kk]=j;
					kk++;
				}
	//GetRNGstate();
	time_t  tin,tfin;
	tin = time (NULL);
	for(n=0;n<max_iter;t++)
	{
		//random rewiring
		if(verbose==1)
    	loadBar( n,  max_iter, 100,  50);
		rand1=(size_t) (unif_rand()*e);
		do  rand2=(size_t) (unif_rand()*e);while (rand1==rand2);
		a=from[rand1];
		c=from[rand2];
		b=to[rand1];
		d=to[rand2];
		if(a!=c && d!=b && matrix[a][d]==0 && matrix[c][b]==0)
			{
				to[rand1]=d;
				to[rand2]=b;
				matrix[a][d]=matrix[c][b]=1;
				matrix[a][b]=matrix[c][d]=0;
				n++;
			}
	if(t>MAXITER)
		{
				//PutRNGstate();
				tfin = time (NULL);
				if(verbose==1)
  		 	 printf("DONE in %d seconds \n",-(tin-tfin));

			return (-1);

		}
	}
	//PutRNGstate();
	tfin = time (NULL);
	if(verbose==1)
    printf("DONE in %d seconds \n",-(tin-tfin));


return 0;

}


size_t inline static check(size_t *pos,size_t *to,size_t *index,size_t a,size_t  b,size_t  c,size_t d,size_t e,size_t nr)
{
    size_t i;

    for(i=index[pos[a]];i<index[pos[a]+1];i++)
    {
        if(to[i]==d)
            return(0);
    }

    for(i=index[pos[c]];i<index[pos[c]+1];i++)
        if(to[i]==b)
            return(0);


    return(1);
}

size_t rewire_sparse_bipartite(size_t *from,size_t *to,size_t nc,size_t nr,size_t max_iter,size_t e,size_t verbose)
{
	size_t i,j,kk,n,rand1,rand2;
 	size_t a,b,c,d;
	size_t *index;
	size_t *pos;
	do	index = (size_t*)calloc( (nr+1),sizeof(size_t ));	while(index==NULL);
	do	pos = (size_t*)calloc(e,sizeof(size_t ));	 while(pos==NULL);
	index[0]=0;
  pos[0]=0;
	kk=1;
  j=0;
	for( i=1;i<e;i++)
    {
			if(from[i]!=from[i-1])
        {
					index[kk++]=i;
          j++;
        }
        pos[i]=j;
    }
    index[nr]=e;
	//GetRNGstate();
	time_t  tin,tfin;
	tin = time (NULL);

	for(n=0;n<max_iter;n++)
	{

		//random rewiring
  	if(verbose==1)
			loadBar( n,  max_iter, 100,  50);
		rand1=(size_t) (unif_rand()*e);
		do rand2=(size_t) (unif_rand()*e); while (rand1==rand2);
		a=from[rand1];
		c=from[rand2];
		b=to[rand1];
		d=to[rand2];
		if(a!=c && d!=b &&  check(pos,to,index,rand1,b,rand2,d,e,nr)==1 )
        {
						to[rand1]=d;
            to[rand2]=b;

            //printf("Rewire ok!\n");
        }

	}
	tfin = time (NULL);
	//PutRNGstate();
	if(verbose==1)
		printf("DONE in %d seconds \n",-(tin-tfin));
	return 0;


}


size_t rewire_sparse_bipartite_ex(size_t *from,size_t *to,size_t nc,size_t nr,size_t max_iter,size_t e,size_t verbose,size_t MAXITER)
{
	size_t i,j,kk,n,rand1,rand2,t=0;
 	size_t a,b,c,d;
	size_t *index;
	size_t *pos;
	do	index = (size_t*)calloc( (nr+1),sizeof(size_t ));	while(index==NULL);
	do	pos = (size_t*)calloc(e,sizeof(size_t ));	 while(pos==NULL);
	index[0]=0;
  pos[0]=0;
	kk=1;
  j=0;
	for( i=1;i<e;i++)
    {
			if(from[i]!=from[i-1])
        {
					index[kk++]=i;
          j++;
        }
        pos[i]=j;
    }
    index[nr]=e;
	//GetRNGstate();
	time_t  tin,tfin;
	tin = time (NULL);
	for(n=0;n<max_iter;t++)
	{
		//random rewiring
  	if(verbose==1)
			loadBar( n,  max_iter, 100,  50);
		rand1=(size_t) (unif_rand()*e);
		do rand2=(size_t) (unif_rand()*e); while (rand1==rand2);
		a=from[rand1];
		c=from[rand2];
		b=to[rand1];
		d=to[rand2];
		if(a!=c && d!=b &&  check(pos,to,index,rand1,b,rand2,d,e,nr)==1 )
        {
						to[rand1]=d;
            to[rand2]=b;
            n++;
            //printf("Rewire ok!\n");
        }
		if(t>MAXITER)
		{
			tfin = time (NULL);
			//PutRNGstate();
			if(verbose==1)
				printf("DONE in %d seconds \n",-(tin-tfin));

			return (-1);
		}
	}
	tfin = time (NULL);
	//PutRNGstate();
	if(verbose==1)
		printf("DONE in %d seconds \n",-(tin-tfin));
	return 0;


}

double similarity_undirected(short **m,short **n,size_t ncol,size_t nrow,size_t e)
{

    size_t num=0;
    size_t i, j;
    for(i=0;i<nrow;++i)
        for(j=0;j<i;++j)
            num += m[i][j] * n[i][j];
    return((double)num/(2.0*e-num));
}



size_t analysis_undirected(short **incidence,size_t ncol, size_t nrow,double **scores,size_t step,size_t max_iter,size_t verbose)
{


	size_t i,j,kk,n,index,rand1,rand2;
	size_t dim=(size_t)round((double)max_iter/step)+2;
	size_t e=0;
    //copy of the original incidence matrix
  short **matrix;
	size_t *from;
	size_t *to;
	size_t a,b,c,d;

	do	matrix=(short **)calloc(nrow,sizeof(short*));	while(matrix==NULL);
	for(i=0;i<nrow;i++)
		{
			do	matrix[i]= (short*)calloc(ncol,sizeof(short));	 while(matrix[i]==NULL);
			for(j=0;j<ncol;j++)
        {
            matrix[i][j]=incidence[i][j];
            e+=incidence[i][j];
        }
	}
	e/=2;
	//initialization of score vector overwriting the original TO CHECK
	do	from=(size_t*)calloc(e,sizeof(size_t));	while(from==NULL);
	do 	to=(size_t*)calloc(e,sizeof(size_t));    while(to==NULL);
	kk=0;
	for(i=0;i<nrow;++i)
		for(j=0;j<i;++j)
			if(matrix[i][j]==1)
				{
                from[kk]=i;
                to[kk]=j;
                kk++;
				}
	do	*scores=(double*)calloc(dim,sizeof(double));  while(scores==NULL);
	for(i=0;i<dim;(*scores)[i++]=0.0);
	(*scores)[0]=1.0;
	//RANDOM GENERATOR (MILLISECONDS)
	//GetRNGstate();
	index=1;
	time_t  tin,tfin;
	tin = time (NULL);
	for(n=0;n<max_iter;n++)
	{
		if(verbose==1)
			loadBar( n,  max_iter, 100,  50);
		//random rewiring

		rand1=(size_t) (unif_rand()*e);
		do	rand2=(size_t) (unif_rand()*e);	while (rand1==rand2);
		a=from[rand1];
		c=from[rand2];
		b=to[rand1];
		d=to[rand2];

		//printf("a=%d b=%d c=%d d=%d\n ",a+1,b+1,c+1,d+1);
		if( a!=d && c!=b &&
           (	(incidence[a][d]==0 && incidence[c][b]==0  ) ||
            (incidence[a][c]==0 && incidence[d][b]==0  ) ))
        {
            if(incidence[a][d]==0 && incidence[c][b]==0 && incidence[a][c]==0 && incidence[d][b]==0 )
            {
                if(unif_rand()>=0.5)
                {
                    incidence[a][d]=1;incidence[d][a]=1;
                    incidence[c][b]=1;incidence[b][c]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=d;
                    to[rand2]=b;

                }
                else
                {
                    incidence[a][c]=1;incidence[c][a]=1;
                    incidence[d][b]=1;incidence[b][d]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    //	from[rand1]=d;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;


                }
            }
            else
                if(incidence[a][d]==0 && incidence[c][b]==0)
                {
                    incidence[a][d]=1;incidence[d][a]=1;
                    incidence[c][b]=1;incidence[b][c]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=d;
                    to[rand2]=b;

                }
                else
                {
                    incidence[a][c]=1;incidence[c][a]=1;
                    incidence[d][b]=1;incidence[b][d]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;


                }
        }

   		if(n%step==0)
            {
								(*scores)[index++]=similarity_undirected(matrix,incidence,ncol,nrow,e);
										//printf("%lf \n",(*scores)[index-1]);
						}

	}

	tfin = time (NULL);
	if(verbose==1)
    printf("DONE in %d seconds \n",-(tin-tfin));

  //PutRNGstate();
	return (index-1);
}

size_t analysis_undirected_ex(short **incidence,size_t ncol, size_t nrow,double **scores,size_t step,size_t max_iter,size_t verbose,size_t MAXITER)
{


	size_t i,j,kk,n,index,rand1,rand2,t=0;
	size_t dim=(size_t)round((double)max_iter/step)+2;
	size_t e=0;
    //copy of the original incidence matrix
  short **matrix;
	size_t *from;
	size_t *to;
	size_t a,b,c,d;

	do	matrix=(short **)calloc(nrow,sizeof(short*));	while(matrix==NULL);
	for(i=0;i<nrow;i++)
		{
			do	matrix[i]= (short*)calloc(ncol,sizeof(short));	 while(matrix[i]==NULL);
			for(j=0;j<ncol;j++)
        {
            matrix[i][j]=incidence[i][j];
            e+=incidence[i][j];
        }
	}
	e/=2;
	//initialization of score vector overwriting the original TO CHECK
	do	from=(size_t*)calloc(e,sizeof(size_t));	while(from==NULL);
	do 	to=(size_t*)calloc(e,sizeof(size_t));    while(to==NULL);
	kk=0;
	for(i=0;i<nrow;++i)
		for(j=0;j<i;++j)
			if(matrix[i][j]==1)
				{
                from[kk]=i;
                to[kk]=j;
                kk++;
				}
	do	*scores=(double*)calloc(dim,sizeof(double));  while(scores==NULL);
	for(i=0;i<dim;(*scores)[i++]=0.0);
	(*scores)[0]=1.0;
	//RANDOM GENERATOR (MILLISECONDS)
	//GetRNGstate();
	index=1;
	time_t  tin,tfin;
	tin = time (NULL);
	for(n=0;n<max_iter;t++)
	{
		if(verbose==1)
			loadBar( n,  max_iter, 100,  50);
		//random rewiring

		rand1=(size_t) (unif_rand()*e);
		do	rand2=(size_t) (unif_rand()*e);	while (rand1==rand2);
		a=from[rand1];
		c=from[rand2];
		b=to[rand1];
		d=to[rand2];
		if(t>MAXITER)
		{
			tfin = time (NULL);
			//PutRNGstate();
			if(verbose==1)
  	  printf("DONE in %d seconds \n",-(tin-tfin));

			printf("Reached the maximum number admissible of iterations!\n");
			return (index-1);
		}
		//printf("a=%d b=%d c=%d d=%d\n ",a+1,b+1,c+1,d+1);
		if( a!=d && c!=b &&
           (	(incidence[a][d]==0 && incidence[c][b]==0  ) ||
            (incidence[a][c]==0 && incidence[d][b]==0  ) ))
        {
            if(incidence[a][d]==0 && incidence[c][b]==0 && incidence[a][c]==0 && incidence[d][b]==0 )
            {
                if(unif_rand()>=0.5)
                {
                    incidence[a][d]=1;incidence[d][a]=1;
                    incidence[c][b]=1;incidence[b][c]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=d;
                    to[rand2]=b;
										n++	; if(n%step==0)
            {
								(*scores)[index++]=similarity_undirected(matrix,incidence,ncol,nrow,e);
										//printf("%lf \n",(*scores)[index-1]);
						}
                }
                else
                {
                    incidence[a][c]=1;incidence[c][a]=1;
                    incidence[d][b]=1;incidence[b][d]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    //	from[rand1]=d;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;
                    n++	; if(n%step==0)
            {
								(*scores)[index++]=similarity_undirected(matrix,incidence,ncol,nrow,e);
										//printf("%lf \n",(*scores)[index-1]);
						}

                }
            }
            else
                if(incidence[a][d]==0 && incidence[c][b]==0)
                {
                    incidence[a][d]=1;incidence[d][a]=1;
                    incidence[c][b]=1;incidence[b][c]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=d;
                    to[rand2]=b;
                    n++	; if(n%step==0)
            {
								(*scores)[index++]=similarity_undirected(matrix,incidence,ncol,nrow,e);
										//printf("%lf \n",(*scores)[index-1]);
						}
                }
                else
                {
                    incidence[a][c]=1;incidence[c][a]=1;
                    incidence[d][b]=1;incidence[b][d]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;
                   n++	;
						if(n%step==0)
            {
								(*scores)[index++]=similarity_undirected(matrix,incidence,ncol,nrow,e);
										//printf("%lf \n",(*scores)[index-1]);
						}
                }
        }



	}

	tfin = time (NULL);
	if(verbose==1)
    printf("DONE in %d seconds \n",-(tin-tfin));

  //PutRNGstate();
	return (index-1);
}


size_t rewire(short **incidence,size_t ncol, size_t nrow,size_t max_iter,size_t verbose)
{

 	size_t i,j,kk,n,rand1,rand2;
	size_t e=0;
    //copy of the original incidence matrix
	size_t *from;
	size_t *to;
	size_t a,b,c,d;


	for(i=0;i<nrow;i++)
		for(j=0;j<ncol;j++)
			e+=incidence[i][j];


	e/=2;
	//initialization of score vector overwriting the original TO CHECK
	do from=(size_t*)calloc(e,sizeof(size_t)); while(from==NULL);
	do to=(size_t*)calloc(e,sizeof(size_t));   while(to==NULL);
	kk=0;
	for(i=0;i<nrow;++i)
		for(j=0;j<i;++j)
			if(incidence[i][j]==1)
				{
					from[kk]=i;
					to[kk]=j;
          kk++;
        }
	time_t  tin,tfin;
  tin = time (NULL);
 	//GetRNGstate();
	for(n=0;n<max_iter;n++)
		{
			if(verbose==1)
				loadBar( n,  max_iter, 100,  50);
		//random rewiring
		  rand1=(size_t) (unif_rand()*e);
      do rand2=(size_t) (unif_rand()*e); while (rand1==rand2);
      a=from[rand1];
      c=from[rand2];
      b=to[rand1];
      d=to[rand2];

		//  printf("%d %d %d %d %d %d %d %d %d %d\n ",rand1,rand2,a+1,b+1,c+1,d+1,incidence[a][d],incidence[c][b],incidence[a][c],incidence[d][b]);
     if(a!=c && b!=d&&  a!=d && c!=b &&
           (	(incidence[a][d]==0 && incidence[c][b]==0  ) ||
            (incidence[a][c]==0 && incidence[d][b]==0  ) ))
        {
            if(incidence[a][d]==0 && incidence[c][b]==0 && incidence[a][c]==0 && incidence[d][b]==0 )
            {
                if(unif_rand()>=0.5)
                {
                    incidence[a][d]=1;incidence[d][a]=1;
                    incidence[c][b]=1;incidence[b][c]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=d;
                    to[rand2]=b;


                }
                else
                {
                    incidence[a][c]=1;incidence[c][a]=1;
                    incidence[d][b]=1;incidence[b][d]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    //	from[rand1]=d;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;



                }
            }
            else
                if(incidence[a][d]==0 && incidence[c][b]==0)
                {
                    incidence[a][d]=1;incidence[d][a]=1;
                    incidence[c][b]=1;incidence[b][c]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=d;
                    to[rand2]=b;

                }
                else
                {
                    incidence[a][c]=1;incidence[c][a]=1;
                    incidence[d][b]=1;incidence[b][d]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;

                }
        }
    }
	tfin = time (NULL);
 	if(verbose==1)
    printf("DONE in %d seconds \n",-(tin-tfin));
	//PutRNGstate();

	return 0;


}


size_t rewire_ex(short **incidence,size_t ncol, size_t nrow,size_t max_iter,size_t verbose,size_t MAXITER)
{

 	size_t i,j,kk,n,rand1,rand2,t=0;
	size_t e=0;
    //copy of the original incidence matrix
	size_t *from;
	size_t *to;
	size_t a,b,c,d;


	for(i=0;i<nrow;i++)
		for(j=0;j<ncol;j++)
			e+=incidence[i][j];


	e/=2;
	//initialization of score vector overwriting the original TO CHECK
	do from=(size_t*)calloc(e,sizeof(size_t)); while(from==NULL);
	do to=(size_t*)calloc(e,sizeof(size_t));   while(to==NULL);
	kk=0;
	for(i=0;i<nrow;++i)
		for(j=0;j<i;++j)
			if(incidence[i][j]==1)
				{
					from[kk]=i;
					to[kk]=j;
          kk++;
        }
	time_t  tin,tfin;
  tin = time (NULL);

 	//GetRNGstate();
	for(n=0;n<max_iter;t++)
		{
			if(verbose==1)
				loadBar( n,  max_iter, 100,  50);
		//random rewiring
		  rand1=(size_t) (unif_rand()*e);
      do rand2=(size_t) (unif_rand()*e); while (rand1==rand2);
      a=from[rand1];
      c=from[rand2];
      b=to[rand1];
      d=to[rand2];
	if(t>MAXITER)
			{
				tfin = time (NULL);
 				if(verbose==1)
    		printf("DONE in %d seconds \n",-(tin-tfin));
				//PutRNGstate();
				return (-1);

			}
		//  printf("%d %d %d %d %d %d %d %d %d %d\n ",rand1,rand2,a+1,b+1,c+1,d+1,incidence[a][d],incidence[c][b],incidence[a][c],incidence[d][b]);
     if(a!=c && b!=d&&  a!=d && c!=b &&
           (	(incidence[a][d]==0 && incidence[c][b]==0  ) ||
            (incidence[a][c]==0 && incidence[d][b]==0  ) ))
        {
            if(incidence[a][d]==0 && incidence[c][b]==0 && incidence[a][c]==0 && incidence[d][b]==0 )
            {
                if(unif_rand()>=0)
                {
                    incidence[a][d]=1;incidence[d][a]=1;
                    incidence[c][b]=1;incidence[b][c]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=d;
                    to[rand2]=b;
                    n++;

                }
                else
                {
                    incidence[a][c]=1;incidence[c][a]=1;
                    incidence[d][b]=1;incidence[b][d]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    //	from[rand1]=d;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;

                    n++;

                }
            }
            else
                if(incidence[a][d]==0 && incidence[c][b]==0)
                {
                    incidence[a][d]=1;incidence[d][a]=1;
                    incidence[c][b]=1;incidence[b][c]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=d;
                    to[rand2]=b;
                    n++;
                }
                else
                {
                    incidence[a][c]=1;incidence[c][a]=1;
                    incidence[d][b]=1;incidence[b][d]=1;
                    incidence[a][b]=0;incidence[b][a]=0;
                    incidence[c][d]=0;incidence[d][c]=0;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;
                    n++;
                }
        }
    }
	tfin = time (NULL);
 	if(verbose==1)
    printf("DONE in %d seconds \n",-(tin-tfin));
	//PutRNGstate();

	return 0;


}

void static inline sub(size_t a,size_t b,size_t c,size_t d,size_t *degree,short **adj)
{
    size_t i;

    for(i=0;i<degree[a];i++)
        if(adj[a][i]==b)
        {
            adj[a][i]=d;
            break;
        }

    for(i=0;i<degree[b];i++)
        if(adj[b][i]==a)
        {
            adj[b][i]=c;
            break;
        }
    for(i=0;i<degree[c];i++)
        if(adj[c][i]==d)
        {
            adj[c][i]=b;
            break;
        }

    for(i=0;i<degree[d];i++)
        if(adj[d][i]==c)
        {
            adj[d][i]=a;
            break;
        }


}

void static inline sub2(size_t a,size_t b,size_t c,size_t d,size_t *degree,short **adj)
{
    size_t i;

    for(i=0;i<degree[a];i++)
        if(adj[a][i]==b)
        {
            adj[a][i]=c;
            break;
        }

    for(i=0;i<degree[b];i++)
        if(adj[b][i]==a)
        {
            adj[b][i]=d;
            break;
        }
    for(i=0;i<degree[c];i++)
        if(adj[c][i]==d)
        {
            adj[c][i]=a;
            break;
        }

    for(i=0;i<degree[d];i++)
        if(adj[d][i]==c)
        {
            adj[d][i]=b;
            break;
        }


}

size_t static inline is_not(size_t a,size_t d,size_t *degree,short **adj)
{
    size_t i;

    for(i=0;i<degree[a];i++)
        if(adj[a][i]==d)
            return(0);
    return(1);
}
size_t rewire_sparse(size_t *from, size_t *to,size_t *degree,size_t ncol, size_t nrow,size_t max_iter, size_t e,size_t verbose)
{

 	size_t i,n,rand1,rand2;
    //copy of the original incidence matrix
	size_t a,b,c,d;
	size_t ad,cb,ac,bd;
	short **adj;
	do	adj=(short **)calloc(nrow,sizeof(short*)); while(adj==NULL);
	for(i=0;i<nrow;i++)
    {
       do	adj[i]= (short*)calloc(degree[i]+1,sizeof(short)); while(adj[i]==NULL);
        adj[i][degree[i]]= degree[i];
    }
	for(i=0;i<e;i++)
		{
			a=from[i];
			b=to[i];
			adj[a][degree[a]-adj[a][degree[a]]]=b;
			adj[a][degree[a]]--;
			adj[b][degree[b]-adj[b][degree[b]]]=a;
			adj[b][degree[b]]--;
    }

  //GetRNGstate();

  time_t  tin,tfin;
	tin = time (NULL);
	for(n=0;n<max_iter;n++)
	{
		if(verbose==1)
			loadBar( n,  max_iter, 100,  50);
		//random rewiring
		   rand1=(size_t) (unif_rand()*e);
        do   rand2=(size_t) (unif_rand()*e); while (rand1==rand2);
        a=from[rand1];
        c=from[rand2];
        b=to[rand1];
        d=to[rand2];

        ad=is_not(a,d,degree,adj);
        cb=is_not(c,b,degree,adj);
        ac=is_not(c,a,degree,adj);
        bd=is_not(b,d,degree,adj);

		  // printf("%d %d %d %d %d %d %d %d %d %d \n ",rand1, rand2,a+1,b+1,c+1,d+1,ad,cb,ac,bd);

        if( a!=c && b!=d & a!=d && c!=b && (( ad==1 && cb==1 )|| (ac==1 && bd==1) ))
        {
            if(ad==1 && cb==1 && ac==1 && bd==1 )
            {
                if(unif_rand()>=0.5)
                {

                    to[rand1]=d;
                    to[rand2]=b;
                    sub(a,b,c,d,degree,adj);


                }
                else
                {
                    sub2(a,b,c,d,degree,adj);

                    //	from[rand1]=d;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;


                }
            }
            else
                if(ad==1 && cb==1 )
                {

                    to[rand1]=d;
                    to[rand2]=b;
                    sub(a,b,c,d,degree,adj);



                }
                else
                {
                    sub2(a,b,c,d,degree,adj);


                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;

                }
        }


    }
  tfin = time (NULL);
  if(verbose==1)
 	 printf("DONE in %d seconds \n",-(tin-tfin));
  //PutRNGstate();
	return 0;
}

size_t rewire_sparse_ex(size_t *from, size_t *to,size_t *degree,size_t ncol, size_t nrow,size_t max_iter, size_t e,size_t verbose,size_t MAXITER)
{

 	size_t i,n,rand1,rand2,t=0;
    //copy of the original incidence matrix
	size_t a,b,c,d;
	size_t ad,cb,ac,bd;
	short **adj;
	do	adj=(short **)calloc(nrow,sizeof(short*)); while(adj==NULL);
	for(i=0;i<nrow;i++)
    {
       do	adj[i]= (short*)calloc(degree[i]+1,sizeof(short)); while(adj[i]==NULL);
        adj[i][degree[i]]= degree[i];
    }
	for(i=0;i<e;i++)
		{
			a=from[i];
			b=to[i];
			adj[a][degree[a]-adj[a][degree[a]]]=b;
			adj[a][degree[a]]--;
			adj[b][degree[b]-adj[b][degree[b]]]=a;
			adj[b][degree[b]]--;
    }

  //GetRNGstate();

  time_t  tin,tfin;
	tin = time (NULL);
	for(n=0;n<max_iter;t++)
	{
		if(verbose==1)
			loadBar( n,  max_iter, 100,  50);
		//random rewiring
		   rand1=(size_t) (unif_rand()*e);
        do   rand2=(size_t) (unif_rand()*e); while (rand1==rand2);
        a=from[rand1];
        c=from[rand2];
        b=to[rand1];
        d=to[rand2];

        ad=is_not(a,d,degree,adj);
        cb=is_not(c,b,degree,adj);
        ac=is_not(c,a,degree,adj);
        bd=is_not(b,d,degree,adj);
			if(t>MAXITER)
		{

  tfin = time (NULL);
	if(verbose==1)
    printf("DONE in %d seconds \n",-(tin-tfin));
   //PutRNGstate();
			return (-1);




			}
		  // printf("%d %d %d %d %d %d %d %d %d %d \n ",rand1, rand2,a+1,b+1,c+1,d+1,ad,cb,ac,bd);

        if( a!=c && b!=d & a!=d && c!=b && (( ad==1 && cb==1 )|| (ac==1 && bd==1) ))
        {
            if(ad==1 && cb==1 && ac==1 && bd==1 )
            {
                if(unif_rand()>=0.5)
                {

                    to[rand1]=d;
                    to[rand2]=b;
                    sub(a,b,c,d,degree,adj);

                    n++;
                }
                else
                {
                    sub2(a,b,c,d,degree,adj);

                    //	from[rand1]=d;
                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;
                               n++;

                }
            }
            else
                if(ad==1 && cb==1 )
                {

                    to[rand1]=d;
                    to[rand2]=b;
                    sub(a,b,c,d,degree,adj);


                    n++;
                }
                else
                {
                    sub2(a,b,c,d,degree,adj);


                    to[rand1]=c;
                    from[rand2]=b;
                    to[rand2]=d;
                    n++;
                }
        }


    }
  tfin = time (NULL);
  if(verbose==1)
 	 printf("DONE in %d seconds \n",-(tin-tfin));
  //PutRNGstate();
	return 0;
}
